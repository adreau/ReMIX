#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
using namespace std;

const char tab = '\t';

vector<int> variantPosGlobal;
vector<int> variantQualGlobal;
vector<string> genoBarcodes_prePhasing;
map<int,string> genoBarcodes_postPhasing;

map<int,int> variantH0MolCoverage;
map<int,int> variantH1MolCoverage;
map<int,int> variantHmixMolCoverage;
map<int,int> variantHotherMolCoverage;

vector<int> variantMolCoverage;

float meanMolCoverage;

string lastLinePrevChrPre="";

/*recompute h0, h1 and hmix probabilities
  parameters: fragment_phasing.tsv
            Populate info fields vcf file
            Phase snp indels vcf file
*/

void passHeader(ifstream &file){

  string line;
  getline(file,line);

  while(line[0]=='#' && line[1]=='#' ){
    getline(file,line);
  }


}


//get the informations for a molecule from a line in fragment_phasing file
void getMoleculeInfo( string& line, string& chr,
                      int& start, int& end, string& barcode, int& noReads){

      istringstream ssRead(line);
      string tmp;

      getline(ssRead, chr, tab);

      getline(ssRead, tmp, tab);
      start = stoi(tmp);

      getline(ssRead, tmp, tab);
      end = stoi(tmp);

      getline(ssRead, barcode, tab);

      getline(ssRead, tmp, tab);
      noReads = stoi(tmp);

}

//get variant information from a file in the vcf file
void getVariantInfo( string &chr, int& vPos, int& vQual,
                    string& genoBC_pre){

    istringstream ssLinePre(lastLinePrevChrPre);

    string tmp;

    //col 1 chr
    getline(ssLinePre, chr, tab);

    //col 2 variant position
    getline(ssLinePre, tmp, tab);
    vPos=stoi(tmp);

    //col 3 to 5: ID, REF, ALT
    for(int i=0;i<3;i++){
      getline(ssLinePre, tmp, tab);
    }

    //col 6 variant quality
    getline(ssLinePre, tmp, tab);
    vQual=stoi(tmp);

    //col 7 8 9: FILTER, INFO and FORMAT
    for(int i=0;i<3;i++){
      getline(ssLinePre, tmp, tab);
    }

    //col 10: genotype and BC info
    getline(ssLinePre, genoBC_pre, tab);

}

void chargeNextChr( string& chrCurent, ifstream &vcfPrePhasing ){

  variantPosGlobal.clear();
  variantQualGlobal.clear();
  genoBarcodes_prePhasing.clear();


  string chr, genoBC_pre;
  int vPos, vQual;

  if(lastLinePrevChrPre.length() == 0){

    getline(vcfPrePhasing, lastLinePrevChrPre);

  }

  //cout<<lastLinePrevChrPre<<endl;
  //cout<<lastLinePrevChrPost<<endl;

  getVariantInfo( chr, vPos, vQual, genoBC_pre);

  while( chrCurent.compare(chr) != 0 && getline(vcfPrePhasing, lastLinePrevChrPre) ){

    getVariantInfo( chr, vPos, vQual, genoBC_pre);

  }

  //cout<<"chargeNextChr: "<<chr<<" "<<vPos<<" "<<vQual<<" "<<genoBC_pre<<" "<<genoBC_post<<endl;

  while (chrCurent.compare(chr) == 0) {

    variantPosGlobal.push_back(vPos);
    variantQualGlobal.push_back(vQual);
    genoBarcodes_prePhasing.push_back(genoBC_pre);


    if(getline(vcfPrePhasing, lastLinePrevChrPre) )
      getVariantInfo(chr, vPos, vQual, genoBC_pre);
    else break;

  }

}

void chargePostPhasingInfo(ifstream &postPhasingFile){

  genoBarcodes_postPhasing.clear();

  //initialise with the prephasing info
  /*for(int i=0; i<genoBarcodes_prePhasing.size(); i++)
    genoBarcodes_postPhasing.insert( pair<int,string>
        (variantPosGlobal[i],genoBarcodes_prePhasing[i].substr(0,3) ) );
*/
  string line;
  int vPos;
  string geno1;
  string geno2;
  string tmp;
  string genoPost;
  map<int,string>::iterator itGenoPost;
  //read HBOP file
  while(getline(postPhasingFile,line)){

    if(line[0]!='B' && line[0]!='*'){
      istringstream ssLinePost(line);
      getline(ssLinePost, tmp, tab);
      vPos = stoi(tmp);


      getline(ssLinePost, geno1, tab);
      getline(ssLinePost, geno2, tab);


      if(geno1.compare("-")!=0 && geno2.compare("-")!=0){

        genoPost=geno1+"|"+geno2;
        //genoBarcodes_postPhasing.insert( pair<int,string>(vPos,genoPost);
        genoBarcodes_postPhasing[vPos]=genoPost;

      }
    }


  }

}

//get all the variants info for a barcode, a start and an end
void getBCInfo(string& barcode, int& start, int& end, vector<int>& varPos,
  vector<double>& var_quality, vector<int>& genoPre, vector<int>& genoPost,
  vector<string>& bcQual_h0, vector<string>& bcQual_h1){


    int varIndex=0;
    string genoStr_pre="";
    string genoStr_post="";

    int geno1, geno2;
    int hap1, hap2;
    bool phased;

    int posBC1, posBC2;
    int posComma;
    int posBC1_end, posBC2_end;

    int prevVarPos=0;
    int maxQual=0;
    int varPosCurr=0;

    string qualsBC1,qualsBC2;

    int variantPosGlobalSize=variantPosGlobal.size();

    while(variantPosGlobal[varIndex]<start && varIndex<variantPosGlobalSize-1){
      varIndex++;
    }

    //cout<<varIndex<<" "<<variantPosGlobal.size()<<endl;

    while(variantPosGlobal[varIndex]<=end && varIndex<variantPosGlobalSize){

      varPosCurr=variantPosGlobal[varIndex];

      genoStr_pre=genoBarcodes_prePhasing[varIndex];

      posBC1 = genoStr_pre.find(barcode);

      //if the barcode span the variant with at least a read and
      //if the variant was not filtered and was phased by HBOP
      if( posBC1 != string::npos &&
        genoBarcodes_postPhasing.find(varPosCurr) != genoBarcodes_postPhasing.end()){

        //cout<<varIndex<<endl;
        //cout<<genoStr_pre<<endl;
        //cout<<genoStr_post<<endl;

        genoStr_post=genoBarcodes_postPhasing.at(varPosCurr);
        varPos.push_back(varPosCurr);
        var_quality.push_back(variantQualGlobal[varIndex]);


        //treat the pre phasing info
        geno1 = 0;
        geno2 = 1;

        posComma = genoStr_pre.find(',');

        posBC1_end=genoStr_pre.find(";",posBC1);
        if(posBC1_end > posComma && posComma > posBC1) posBC1_end=posComma;
        qualsBC1=genoStr_pre.substr(posBC1+19,posBC1_end-posBC1-19);


        posBC2 = genoStr_pre.find(barcode,posBC1+1);

        //cout<<posBC1<<" " <<posBC2<<endl;

        if(posBC2 == string::npos){

          qualsBC2="";

          if(posBC1<posComma){

            genoPre.push_back(geno1);
            bcQual_h0.push_back(qualsBC1);
            bcQual_h1.push_back(qualsBC2);

          }else{

            genoPre.push_back(geno2);
            bcQual_h0.push_back(qualsBC2);
            bcQual_h1.push_back(qualsBC1);

          }


        }else{
            genoPre.push_back(2);

            posBC2_end=genoStr_pre.find(";",posBC2);
            qualsBC2=genoStr_pre.substr(posBC2+19,posBC2_end-posBC2-19);

            if(geno1 == 0){

                bcQual_h0.push_back(qualsBC1);
                bcQual_h1.push_back(qualsBC2);

            }else{

              bcQual_h0.push_back(qualsBC2);
              bcQual_h1.push_back(qualsBC1);

            }
        }

        //treat the post phasing info
        hap1 = genoStr_post[0] - '0';
        hap2 = genoStr_post[2] - '0';

        if(genoStr_post[1] == '|'){


          if(posBC2==string::npos){

            if(posBC1<posComma){ //spanes the reference
              if(hap1 == 0) //reference is in hap1
                genoPost.push_back(0);
              else //reference is in hap2
                genoPost.push_back(1);
            }else{ //spanes the allele
              if(hap1 == 1) //the alle is in hap1
                genoPost.push_back(0);
              else //reference is in hap2
                genoPost.push_back(1);
            }

          }else{
              genoPost.push_back(2);
          }

        }else{
          genoPost.push_back(-1);
        }


      }

      varIndex++;

    }

}

//compute h0,h1 and hmix for the current barcode, start and end
void computeProbabilities(vector<int>& varPos, vector<double>& var_quality,
  vector<int>& genotype, vector<int>& genotype2,
  vector<string>& BC_quality1, vector<string>& BC_quality2,
  double& h0, double& h1, double& hmix){

  vector<double> hap1_matrix(varPos.size(), 0);
  vector<double> hap2_matrix(varPos.size(), 0);
  vector<double> mix_matrix(varPos.size(), 0);

  vector<int> bc_quals1, bc_quals2;
  vector<double> quals1, quals2;
  //compute probability for each variant

  for (int i=0;i<var_quality.size();i++){

      double qual=var_quality[i];

      qual = min(50.0, max(1.0, qual));
      double prob_var_wrong = pow(10.0, -qual/10.0);

      double log_prob_wrong = log(prob_var_wrong);
      double log_prob_notwrong = log(1.0 - prob_var_wrong);

      int ILLUMINA_QUAL_OFFSET=33;

      istringstream bc_quals1Str(BC_quality1[i]);
      istringstream bc_quals2Str(BC_quality2[i]);

      bc_quals1.clear();
      bc_quals2.clear();

      string bc_qual;
      int bc_qual_i;
      while(getline(bc_quals1Str,bc_qual,'_')){

        //!!filter BC quality
        bc_qual_i=stoi(bc_qual);
        //if(bc_qual_i >= 70)
          bc_quals1.push_back(bc_qual_i);

      }

      while(getline(bc_quals2Str,bc_qual,'_')){

        bc_qual_i=stoi(bc_qual);
        //if(bc_qual_i >= 70)
          bc_quals2.push_back(bc_qual_i);

      }

      //ignore the variant if the barcode has no read with a quality > 70
      /*if(bc_quals1.size()==0 && bc_quals2.size()==0){

        genotype2[i]=-2;

      }*/

      //if we lost all reads covering one allele where there should be 2
      /*if((genotype[i]==2) && ( bc_quals1.size()==0 || bc_quals2.size()==0)){

        genotype[i]=(bc_quals1.size()>0)?0:1;
        genotype2[i]=0;

      }*/

      quals1.clear();
      quals2.clear();

      //if(genotype[i]==0){
      for(int j=0;j<bc_quals1.size();j++)
        quals1.push_back(log(pow(10.0, -(bc_quals1[j] - ILLUMINA_QUAL_OFFSET)/10.0)));
      for(int j=0;j<bc_quals2.size();j++)
        quals2.push_back(log(pow(10.0, -(bc_quals2[j] - ILLUMINA_QUAL_OFFSET)/10.0)));



      int total_length = quals1.size() + quals2.size();

      double sumQuals1=0;
      double logQuals1=0;
      for(int j=0;j<quals1.size();j++){

          sumQuals1+=quals1[j];
          logQuals1+=log(1.0 - exp(quals1[j]));

      }

      double sumQuals2=0;
      double logQuals2=0;
      for(int j=0;j<quals2.size();j++){

          sumQuals2+=quals2[j];
          logQuals2+=log(1.0 - exp(quals2[j]));

      }

      double log_hap_prob_1 = logQuals1 + sumQuals2;
      double log_hap_prob_2 = logQuals2 + sumQuals1;
      double log_hap_prob_mix=0;

      for(int j=0;j<total_length;j++){
        log_hap_prob_mix +=log(0.5);
      }



      hap1_matrix[i] = log( exp(log_prob_wrong) + exp(log_prob_notwrong + log_hap_prob_1));
      hap2_matrix[i] = log( exp(log_prob_wrong) + exp(log_prob_notwrong + log_hap_prob_2));
      mix_matrix[i] = log( exp(log_prob_wrong) + exp(log_prob_notwrong + log_hap_prob_mix));

  }

  int compPos=0;
  int varSize=varPos.size();

  int prevPos=0;
  int maxQual=0;
  int maxQualPos=-1;
  int firstPosInt=0;

  while(compPos < varSize){

    if(varPos[compPos] == prevPos){

      if(maxQual < var_quality[compPos]){

          maxQual = var_quality[compPos];
          maxQualPos = compPos;
      }

    }else{

      if(maxQualPos>=0){

        for(int i=firstPosInt; i<compPos ; i++)
          if(i!=maxQualPos)
            genotype2[i]=-2;

      }
      firstPosInt = compPos;
      maxQualPos = -1;
      maxQual = 0;
    }

    prevPos=varPos[compPos];

    compPos++;

  }

  if(maxQualPos>=0){

    for(int i=firstPosInt; i<compPos ; i++)
      if(i!=maxQualPos)
        genotype2[i]=-2;

  }


  //compute probability for each barcode

  double bc_mix_prob=0.001L;

  double p1 = 0;
  double p2 = 0;
  double p_mix = 0;

  double l1m = log(0.5 * (1.0 - bc_mix_prob));
  double l_mix = log(bc_mix_prob);

  for(int i=0; i<genotype.size();i++) {

    /*genotype2 =-2 => ignore the variant
      genotype2 = -1 => variant not phased in PHASED_VARIANTS.vcf
      genotype2= genotype phasing didn't change
      genotype2 != genotype phasing changed
    */
    if(genotype2[i] != -2){
      if (genotype2[i] == 0){
        if(genotype[i]==0){
          p1 += hap1_matrix[i];
          p2 += hap2_matrix[i];
        }else{
          p1 += hap2_matrix[i];
          p2 += hap1_matrix[i];
        }
      }else if (genotype2[i] == 1) {
        if(genotype[i]==0){
          p1 += hap2_matrix[i];
          p2 += hap1_matrix[i];
        }else{
          p1 += hap1_matrix[i];
          p2 += hap2_matrix[i];
        }
      }else {
          p1 += mix_matrix[i];
          p2 += mix_matrix[i];
      }
      p_mix += mix_matrix[i];
    }
  }

  p1 += l1m;
  p2 += l1m;
  p_mix += l_mix;



  double norm = log( exp(p1) + exp(p2) +exp(p_mix));

  h0 = exp(p1 - norm);
  h1 = exp(p2 - norm);
  hmix = exp(p_mix - norm);


}

float countPhase(vector<int>& genoPost, int posBegin, int posEnd, int phase){

  float sum=0;
  for(int i=posBegin;i<posEnd;i++){
    if(genoPost[i]==phase)
      sum++;
  }


  return sum/(posEnd-posBegin);

}

bool checkLengthOfSwitches(vector<int>& genoPost, vector<int>& variantPos, int& var1, int& var2){

  bool isRealHmixMol=true;
  float per0, per1=0;

  for(int i=1; i< genoPost.size();i++){

    if((genoPost[i]==0)&&(genoPost[i-1]==1)){

        per0=countPhase(genoPost,i,genoPost.size(),0);
        per1=countPhase(genoPost,0,i,1);

    }else if ((genoPost[i]==1)&&(genoPost[i-1]==0)){

        per1=countPhase(genoPost,i,genoPost.size(),1);
        per0=countPhase(genoPost,0,i,0);

    }

    if((per0>=0.7)&&(per1>=0.7)){
      var1=variantPos[i-1];
      var2=variantPos[i];
      return true;
    }

  }

  return false;
}

int main (int argc, char* argv[])
{

  cout.precision(13);

  ifstream fragPh(argv[1]);

  ifstream vcfPrePhasing(argv[2]);

  ifstream postPhasingFile(argv[3]);

  string chromosome = argv[4];

  passHeader(vcfPrePhasing);

  string line="";
  string currentChr="";

  string barcode="";
  string chr="";
  int start;
  int end;

  vector<int> variantPos;
  vector<double> variantQual;
  vector<int> genoPre;
  vector<int> genoPost;
  vector<string> bcQual_h0;
  vector<string> bcQual_h1;

  double h0, h1, hmix;

  int noReads;

  getline(fragPh,line);
  getMoleculeInfo(line, chr, start, end, barcode, noReads);

  while(chr.compare(chromosome) !=0 && getline(fragPh,line)){

      getMoleculeInfo(line, chr, start, end, barcode, noReads);

  }


  while (chr.compare(chromosome) == 0){

    if (currentChr.compare(chr)!=0){

        chargeNextChr(chr,vcfPrePhasing);
        chargePostPhasingInfo(postPhasingFile);

        currentChr=chr;

        for(int i=0; i<variantPosGlobal.size();i++){
          variantH0MolCoverage.insert(pair<int,int>(variantPosGlobal[i],0));
          variantH1MolCoverage.insert(pair<int,int>(variantPosGlobal[i],0));
          variantHmixMolCoverage.insert(pair<int,int>(variantPosGlobal[i],0));
          variantHotherMolCoverage.insert(pair<int,int>(variantPosGlobal[i],0));
        }

    }

    variantPos.clear();
    variantQual.clear();
    genoPre.clear();
    genoPost.clear();
    bcQual_h0.clear();
    bcQual_h1.clear();

    getBCInfo(barcode, start, end, variantPos, variantQual, genoPre, genoPost, bcQual_h0, bcQual_h1);

    if(genoPost.size()>0){
      computeProbabilities(variantPos, variantQual, genoPre, genoPost, bcQual_h0, bcQual_h1, h0, h1, hmix);

      if(h0>h1 && h0>hmix){
        for(int i=0; i<variantPos.size(); i++){
          variantH0MolCoverage.at(variantPos[i])++;
        }
      }
      else if(h1>h0 && h1>hmix){
            for(int i=0; i<variantPos.size(); i++){
              variantH1MolCoverage.at(variantPos[i])++;
            }
          }
      else if(hmix>h0 && hmix>h1){
            for(int i=0; i<variantPos.size(); i++){
              variantHmixMolCoverage.at(variantPos[i])++;
            }
          }
      else
        for(int i=0; i<variantPos.size(); i++){
            variantHotherMolCoverage.at(variantPos[i])++;
        }

    }


    if( getline(fragPh,line) )
        getMoleculeInfo(line, chr, start, end, barcode, noReads);
    else break;


  }

//print haplotype coverage and remove variant if too different
  string varDataName=chromosome;
  varDataName.append(".varMolCov");
  ofstream varData(varDataName.c_str());


  int countVar=0;
  int varPos, h0Mol,h1Mol,hmixMol,hotherMol;
  while( countVar < variantPosGlobal.size() ){

    varPos=variantPosGlobal[countVar];

    h0Mol=variantH0MolCoverage.at(varPos);
    h1Mol=variantH1MolCoverage.at(varPos);
    hmixMol=variantHmixMolCoverage.at(varPos);
    hotherMol=variantHotherMolCoverage.at(varPos);

    float diff = (float)min(h0Mol,h1Mol)/max(h0Mol,h1Mol);
    if(diff<0.82 && hmixMol >= (h0Mol+h1Mol)/10){ //remove the variant

      while(hmixMol >= (h0Mol+h1Mol)/10 && countVar < variantPosGlobal.size()+1){
        variantPosGlobal.erase(variantPosGlobal.begin()+countVar);

        variantQualGlobal.erase(variantQualGlobal.begin()+countVar);
        genoBarcodes_prePhasing.erase(genoBarcodes_prePhasing.begin()+countVar);

        genoBarcodes_postPhasing.erase(varPos);

        varData<<varPos<<tab<<h0Mol<<tab<<h1Mol<<tab<<
              hmixMol<< tab<< hotherMol<<tab<<"Removed"<<endl;

        varPos=variantPosGlobal[countVar];

        h0Mol=variantH0MolCoverage.at(varPos);
        h1Mol=variantH1MolCoverage.at(varPos);
        hmixMol=variantHmixMolCoverage.at(varPos);
        hotherMol=variantHotherMolCoverage.at(varPos);

      }

    }else{

      countVar++;
      varData<<varPos<<tab<<h0Mol<<tab<<h1Mol<<tab<<
            hmixMol<< tab<< hotherMol<<tab<<"Kept"<<endl;

    }

  }

  varData.close();

//re-read fragments and recompute the haplotype probability
  fragPh.clear();
  fragPh.seekg(0, ios::beg);

  getline(fragPh,line);
  getMoleculeInfo(line, chr, start, end, barcode, noReads);

  while(chr.compare(chromosome) !=0 && getline(fragPh,line)){

      getMoleculeInfo(line, chr, start, end, barcode, noReads);

  }

  while (chr.compare(chromosome) == 0){


    variantPos.clear();
    variantQual.clear();
    genoPre.clear();
    genoPost.clear();
    bcQual_h0.clear();
    bcQual_h1.clear();



    getBCInfo(barcode, start, end, variantPos, variantQual, genoPre, genoPost, bcQual_h0, bcQual_h1);

    if(genoPost.size()>0){
        computeProbabilities(variantPos, variantQual, genoPre, genoPost, bcQual_h0, bcQual_h1, h0, h1, hmix);

        if (hmix>h0 && hmix>h1){

          bool xo=0;
          bool mult_xo=0;
          int var1,var2,pos;

          if(checkLengthOfSwitches(genoPost, variantPos, var1, var2)){

            float midPoint=(var1+var2)/2;
            cout<<chr<<tab<<start<<tab<<var1<<tab<<midPoint<<tab<<var2<<tab
                          <<end<<tab<<barcode<<tab<<h0<<tab<<h1<<tab<<hmix<<tab<<noReads<<endl;
          }

      }
    }

    if( getline(fragPh,line) )

        getMoleculeInfo(line, chr, start, end, barcode, noReads);
    else break;
  }

  return 0;

}
