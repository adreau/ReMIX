#!/bin/bash -e

if [ ! -f .envrc ]; then;
    echo 'You need to create a `.envrc` file. See `.envrc.example`.'
    exit 1
fi

source .envrc

testName="$1"
input_R1="$2"
input_R2="$3"
refPath="$4" #10x ref path
refName="$5" #chr3000_1_200000

remix_filter_reads_bin="$remix_binaries_dir"/filterFastq_by_bc
remix_construct_hbop_files_bin="$remix_binaries_dir"/constructHbopFiles_harMean
remix_phase_molecules_bin="$remix_binaries_dir"/XO_fromHBOP_withh0h1parity

mkdir -p ${testName}
cd ${testName}


#
# Step 1 FilterReads
#

"$cutadapt_bin" -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --untrimmed-output ${testName}_R1_001.adapterfiltered.fastq.gz \
    --untrimmed-paired-output ${testName}_R2_001.adapterfiltered.fastq.gz \
    -o ${testName}_R1_001.trimmedadapter.fastq.gz \
    -p ${testName}_R2_001.trimmedadapter.fastq.gz \
    --pair-filter both \
    --info-file ${testName}.R1R2.cutadaptfiltering.info.file \
    ${input_R1} \
    ${input_R2} \
    > output_cutadapt_${testName}.txt 2> error_cutadapt_${testName}.txt

java -jar "$trimmomatic_jar" PE -threads 6 -phred33 \
    ${testName}_R1_001.adapterfiltered.fastq.gz \
    ${testName}_R2_001.adapterfiltered.fastq.gz \
    ${testName}_R1_001.adapterfiltered.avgqual20filtered.fastq.gz ${testName}_R1_001.adapterfiltered.avgqual20filtered.unpaired.fastq.gz \
    ${testName}_R2_001.adapterfiltered.avgqual20filtered.fastq.gz ${testName}_R2_001.adapterfiltered.avgqual20filtered.unpaired.fastq.gz AVGQUAL:20 \
    > output_trimmomatic_${testName}.txt 2> error_trimmomatic_${testName}.txt

mkdir ${testName}
path=$(pwd)
ln -s ${path}/${testName}_R1_001.adapterfiltered.avgqual20filtered.fastq.gz ${testName}/${testName}_S1_L001_R1_001.fastq.gz
ln -s ${path}/${testName}_R2_001.adapterfiltered.avgqual20filtered.fastq.gz ${testName}/${testName}_S1_L001_R2_001.fastq.gz
"$longranger_remix_bin" basic --id=${testName}_basicLongranger --fastqs=${path} --sample=${testName} --localmem=100 --localcores=10 \
    > output_longrangerBasic_${testName}.txt 2> error_longrangerBasic_${testName}.txt

cp ${testName}_basicLongranger/outs/barcoded.fastq.gz ${testName}_R1R2_interleaved.fastq.gz
gunzip ${testName}_R1R2_interleaved.fastq.gz
"$remix_filter_reads_bin" ${testName}_R1R2_interleaved.fastq ${testName}_R1R2_interleaved_filtered_id.fastq \
    > output_filteredReads_${testName}.txt 2> error_filteredReads_${testName}.txt


#
# Step 2 Align reads
#

"$bwa_bin" mem -p -C -M -t 20 \
    "${refPath}/fasta/genome.fa" \
    "${testName}_R1R2_interleaved_filtered_id.fastq" \
    -R "@RG\tID:${testName}\tPL:ILLUMINA\tPU:${testName}.LRsim\tLB:10XGenome\tSM:${testName}" | samtools view -bh - > ${testName}_R1R2.bam

"$samtools_bin" sort -o ${testName}_R1R2_sorted.bam -@ 20 ${testName}_R1R2.bam
"$samtools_bin" index ${testName}_R1R2_sorted.bam

java -Djava.io.tmpdir="$tmp_dir" -Xmx4g -jar "$picard_jar" MarkDuplicates \
    I=${testName}_R1R2_sorted.bam  \
    O=${testName}_R1R2_sorted_mkdup.bam \
    M=${testName}_R1R2_sorted_mkdup_metrics.txt \
    VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE > output_mkdup_${testName}.txt 2> error_mkdup_${testName}.txt

java -Djava.io.tmpdir="$tmp_dir" -Xmx4g -jar "$gatk_jar" \
    -T RealignerTargetCreator \
    -nt 20 \
    -R ${refPath}/fasta/genome.fa \
    -I ${testName}_R1R2_sorted_mkdup.bam \
    -o ${testName}_R1R2_sorted_mkdup.intervals > output_realignTC_${testName}.txt 2> error_realignTC_${testName}.txt

java -Djava.io.tmpdir="$tmp_dir" -Xmx4g -jar "$gatk_jar" \
    -T IndelRealigner \
    -R ${refPath}/fasta/genome.fa \
    -I ${testName}_R1R2_sorted_mkdup.bam \
    -targetIntervals ${testName}_R1R2_sorted_mkdup.intervals \
    -o ${testName}_R1R2_sorted_mkdup_realig.bam > output_indRealign_${testName}.txt 2> error_indRealign_${testName}.txt


#
# Step 3 Variant call
#

echo "${testName}_R1R2_sorted_mkdup_realig.bam" > ${testName}_bam.list

"$samtools_bin" mpileup -b ${testName}_bam.list \
    -A -q 10 -Q 20 -v -u -g \
    --output-tags AD,ADF,ADR,DP,DV,DPR,DP4,SP,INFO/AD,INFO/ADF,INFO/ADR,INFO/DPR \
    -f ${refPath}/fasta/genome.fa | bcftools call -O z -f GQ -M -v -m \
    --threads 20 \
    -o ${testName}.vcf.gz > output_samMPI_${testName}.txt 2> error_samMPI_${testName}.txt

gunzip ${testName}.vcf.gz

"$bcftools_bin" filter -e '(MQSB<0.8)||(MQB<0.4)||(BQB<0.4)||(RPB<0.4)||(FORMAT/GQ<30)||(QUAL<100)||(DP>220)||(DP<5)||(IMF<0.1)||(IMF>0.9)' ${testName}.vcf > ${testName}_filtered.vcf


#
# Step 4 Report molecules
#

"$longranger_remix_bin" mkvcf --reference=${refPath} --sample=${testName} ${testName}_filtered.vcf > output_mkvcf_${testName}.txt 2> error_mkvcf_${testName}.txt

"$samtools_bin" index ${testName}_R1R2_sorted_mkdup_realig.bam

"$samtools_bin" sort -o ${testName}_R1R2_sorted_mkdup_realig_bcSorted.bam -t BX -@ 20 ${testName}_R1R2_sorted_mkdup_realig.bam > \
output_sort_bcSorted_${testName}.txt 2> error_sort_bcSorted_${testName}.txt

"$longranger_remix_bin" reporterMol --id=${testName}_reporterMol --fastqs=test/path --reference=${refPath} \
    --precalled=${path}/sorted_precalled_${testName}.vcf \
    --bam_bc=${path}/${testName}_R1R2_sorted_mkdup_realig_bcSorted.bam \
    --bam_pos=${path}/${testName}_R1R2_sorted_mkdup_realig.bam > output_reporterMol_${testName}.txt 2> error_reporterMol_${testName}.txt

gzip -d -k ${testName}_reporterMol/outs/variants.vcf.gz


#
# Step 5 Phase molecules
#

"$remix_construct_hbop_files_bin" ${testName}_reporterMol/outs/fragments_tsv.tsv \
    ${testName}_reporterMol/outs/variants.vcf ${testName}.frags ${testName}.allVars ${testName}.varData "${refName}" \
    > output_constructHbopFiles_${testName}.txt 2> error_constructHbopFiles_${testName}.txt

(
  cd "$hbop_dir"
  java -cp SIH.jar mpg.molgen.sih.main.SIH -v ${path}/${testName}.allVars -c 2 ${path}/${testName}.frags ${path}/${testName}.phase \
      > ${path}/output_HBOP_${testName}.txt 2> ${path}/error_HBOP_${testName}.txt
)

"$remix_phase_molecules_bin" ${testName}_reporterMol/outs/fragments_tsv.tsv \
    ${testName}_reporterMol/outs/variants.vcf ${testName}.phase "${refName}" > fragment_phasing_vA_HBOP_${testName}.tsv
