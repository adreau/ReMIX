# ReMIX

> Our pipeline called ReMIX is the first method able to identify the recombination rate at a genome wide scale using linked-reads.
> The linked-read information is exploited by ReMIX during three steps: identifying high-quality heterozygous variants, reconstructing molecules, and haplotype phasing each molecule. The molecules identified as recombinant are then used to build an individualized genomic map of recombination cross-overs, enabling us to quantify recombination variation across the genome.


## Dependencies

- g++ >= 4.5
- [longranger-remix](https://github.com/adreau/longranger-remix)
- [Cutadapt](https://github.com/marcelm/cutadapt)
- [Trimmomatic](https://github.com/timflutre/trimmomatic)
- [BWA](https://github.com/lh3/bwa)
- [Picard](https://github.com/broadinstitute/picard)
- [GATK](https://github.com/broadinstitute/gatk)
- [samtools](https://github.com/samtools/samtools)
- [bcftools](https://github.com/samtools/bcftools)
- HBOP : available upon request of the authors of the article [A fast and accurate algorithm for single individual haplotyping](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3521186/)

## Usage

- First configure the paths to the other tools.
```sh
cp .envrc.example .envrc
vim .envrc
```

- Build all the binaries.
```sh
make
```

- Run the pipeline.
```sh
./ReMIX_pipeline.sh
```


## License

The code is licensed under the [GNU Affero General Public License version 3](./LICENSE.md).
