#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Molecule.h"
using namespace std;

const char tab = '\t';

vector<int> variantPosGlobal;
vector<int> variantQualGlobal;
vector<string> genoBarcodes_prePhasing;

vector<Molecule> allMolecules;


vector<int> variantMolCoverage;
vector<int> variantMolWithReadCoverage;

vector<pair<int,int>> variantHapCoverage;

float meanMolCoverage;
float meanMolWithReadCoverage;

string lastLinePrevChrPre="";

/*create H-BOP input files per chromosome
  parameters: fragment_phasing.tsv
            Populate info fields vcf file
            chromosome
*/

//skip vcf header
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

/*consider a Variant only if it's covered by approx the same number of
h0 and h1 molecules that cover the variant with a read*/
bool equivCoverage(string geno_pre){

  bool equiv = true;

  string::iterator it;
  for (it=geno_pre.begin(); it!=geno_pre.end(); ++it)
   if(*it==',')
     break;

  int noMolRef = count (geno_pre.begin(), it, ';');
  noMolRef++;

  int noMolAlt = count (it, geno_pre.end(), ';');
  noMolAlt++;

  if (abs(noMolRef-noMolAlt)>min(noMolRef,noMolAlt)/6)
      equiv=false;
  else
    variantHapCoverage.push_back(pair<int,int>(noMolRef,noMolAlt));

  return equiv;


}

//counts the number of molecules for each variant and computes the mean molecule covereage
void computeMolCovPerVar(){

  int noMolWithRead;

  float totalWithReadMol = 0;
  float totalMol = 0;

  for(int i = 0; i < genoBarcodes_prePhasing.size() ; i++ ){
      noMolWithRead = count (genoBarcodes_prePhasing[i].begin(), genoBarcodes_prePhasing[i].end(), ';');
      noMolWithRead+=2;
      variantMolWithReadCoverage.push_back(noMolWithRead);
      totalWithReadMol+=(float)1/(float)noMolWithRead;

      totalMol+=(float)1/(float)variantMolCoverage[i];

  }

  meanMolWithReadCoverage = (float)genoBarcodes_prePhasing.size()/(float)totalWithReadMol;
  meanMolCoverage = (float)genoBarcodes_prePhasing.size()/(float)totalMol;

}

void removeOvercovVariants(){

  int countVar=0;

  while( countVar < variantPosGlobal.size() ){

    if( ((variantMolCoverage[countVar] > meanMolCoverage*1.5) ||
        (variantMolCoverage[countVar] < meanMolCoverage*0.5)) ||
      ((variantMolWithReadCoverage[countVar] > meanMolWithReadCoverage*1.5)||
      (variantMolWithReadCoverage[countVar] < meanMolWithReadCoverage*0.5)) ){ //remove the variant

      variantPosGlobal.erase(variantPosGlobal.begin()+countVar);

      variantQualGlobal.erase(variantQualGlobal.begin()+countVar);
      genoBarcodes_prePhasing.erase(genoBarcodes_prePhasing.begin()+countVar);
      variantMolCoverage.erase(variantMolCoverage.begin()+countVar);

      variantMolWithReadCoverage.erase(variantMolWithReadCoverage.begin()+countVar);
      variantHapCoverage.erase(variantHapCoverage.begin()+countVar);

    }else{

      countVar++;

    }

  }

}


void chargeNextChrVar( string& chrCurent, ifstream &vcfFile){

  variantPosGlobal.clear();
  variantQualGlobal.clear();
  genoBarcodes_prePhasing.clear();


  string chr, genoBC_pre;
  int vPos, vQual;

  if(lastLinePrevChrPre.length() == 0){

    getline(vcfFile, lastLinePrevChrPre);

  }

  getVariantInfo( chr, vPos, vQual, genoBC_pre);

  while( chrCurent.compare(chr) != 0 && getline(vcfFile, lastLinePrevChrPre)){

    getVariantInfo( chr, vPos, vQual, genoBC_pre);

  }

  while (chrCurent.compare(chr) == 0) {

    if(equivCoverage(genoBC_pre)){
      variantPosGlobal.push_back(vPos);
      variantQualGlobal.push_back(vQual);
      genoBarcodes_prePhasing.push_back(genoBC_pre);

    }
    if(getline(vcfFile, lastLinePrevChrPre))
      getVariantInfo(chr, vPos, vQual, genoBC_pre);
    else break;

    //if(variantPosGlobal.size()%10000==0)
      //cout<<"No of variants charged="<<variantPosGlobal.size()<<endl;

  }

}

void addNewMolToVarCoverage (int start, int end){

  int varIndex=0;
  int variantPosGlobalSize = variantPosGlobal.size();

  while(variantPosGlobal[varIndex]<start && varIndex<variantPosGlobalSize){
    varIndex++;
  }

  //enter in the interval of variants where the pontential spanned variants are
  while(variantPosGlobal[varIndex]<=end && varIndex<variantPosGlobalSize){

    variantMolCoverage[varIndex]++;
    varIndex++;

  }

}


void chargeNextChrFrag( string& chrCurent, ifstream& fragPh){

  //initialise the variant molecule coverage vector
  for(int i=0; i<variantPosGlobal.size();i++){
    variantMolCoverage.push_back(0);
  }

  allMolecules.clear();

  string barcode, chr, line;
  int start,end, noReads;
  getline(fragPh,line);
  getMoleculeInfo(line, chr, start, end, barcode, noReads);

  //go to the right chr in the fragment_phasing file
  while(chr.compare(chrCurent) !=0 && getline(fragPh,line)){

      getMoleculeInfo(line, chr, start, end, barcode, noReads);

  }

  while (chr.compare(chrCurent) == 0){

    if( getline(fragPh,line) ){

      getMoleculeInfo(line, chr, start, end, barcode, noReads);
      Molecule mol(chr, start, end, barcode, noReads);
      allMolecules.push_back(mol);

      //if(allMolecules.size()%10000==0)
        //cout<<"No of molecules charged="<<allMolecules.size()<<endl;

      addNewMolToVarCoverage(start,end);

    }else break;

  }


}

//print info for the molecule
void printMoleculeInfo(Molecule mol, ofstream& fragments){

    string barcode = mol.getBC();
    int start = mol.getStart();
    int end = mol.getEnd();

    int varIndex=0;
    string genoStr_pre="";

    int hap1, hap2;

    int posBC1, posBC2;
    int posComma;

    string lineMolecule;

    int variantPosGlobalSize=variantPosGlobal.size();

    //skip to the variants that can not be spanned by the molecule
    while(variantPosGlobal[varIndex]<start && varIndex<variantPosGlobalSize){
      varIndex++;
    }

    //for each variant of the interval spanned by the molecule
    vector<int> startPositions;
    vector<string> calls;

    int startSegment=-1;
    string callSegment="";

    //enter in the interval of variants where the pontential spanned variants are
    while(variantPosGlobal[varIndex]<=end && varIndex<variantPosGlobalSize){

      //get the barcodes for the current variant
      genoStr_pre=genoBarcodes_prePhasing[varIndex];

      posBC1 = genoStr_pre.find(barcode);

      if( posBC1 != string::npos){//if the molecule span the variant

        if(startSegment == -1)
            startSegment = varIndex+1;
        //  startSegment = variantPosGlobal[varIndex];

        //treat the pre phasing info
        hap1 = genoStr_pre[0] - '0';
        hap2 = genoStr_pre[2] - '0';



        posComma = genoStr_pre.find(',');

        posBC2 = genoStr_pre.find(barcode,posBC1+1);
        //cout<<"Info: "<<hap1<<" "<<hap2<<" "<<posComma<<" "<<posBC1<<endl;
        //cout<<"GenoString: "<<genoStr_pre<<endl;

        //if there is only one allele spanned
        if(posBC2 == string::npos){

          if(posBC1<posComma){
            if(hap1 == 0) //reference is in hap1
              callSegment.append("0");
            else //reference is in hap2
              callSegment.append("1");
          }else{ //spanes the allele
            if(hap1 == 1) //the alle is in hap1
              callSegment.append("0");
            else //reference is in hap2
              callSegment.append("1");
          }
          //cout<<"callSegment: "<<callSegment<<endl;

        }else{ //the molecule has two reads spanning both aleles
          //we cut the molecule in two segments
          if(callSegment.length()>0){ //a segment already started
            startPositions.push_back(startSegment);
            calls.push_back(callSegment);
            callSegment="";
          }
          startSegment=-1;
        }

      }else{//if the molecule does not span the variant

        if(callSegment.length()>0){ //a segment already started

          startPositions.push_back(startSegment);
          calls.push_back(callSegment);
          startSegment=-1;
          callSegment="";

        }

      }

      varIndex++;

    }

    //if the molecule span the last variant of the interval
    if(callSegment.length()>0){

      startPositions.push_back(startSegment);
      calls.push_back(callSegment);
      startSegment=-1;
      callSegment="";

    }


    //print molecules that span at least two fragments
    //if( (startPositions.size()>2)) || (startPositions.size()==1 && calls[0].length() >1 ) ){
    if(startPositions.size()>0){
        //print the line of the molecule
        string moleculeName=barcode + "_" + to_string(start) + "_" +to_string(end);
        lineMolecule="";

        lineMolecule.append(to_string(startPositions.size()));
        lineMolecule.append(string(1,tab));

        lineMolecule.append(moleculeName);


        for(int i = 0; i < startPositions.size(); i++){

          lineMolecule.append(string(1,tab));
          lineMolecule.append(to_string(startPositions[i]));
          lineMolecule.append(string(1,tab));
          lineMolecule.append(calls[i]);

        }

        fragments<<lineMolecule<<endl;
    }


}


void printVarsFile(string chromosome, ofstream& allVars){

  for(int i = 0; i < variantPosGlobal.size() ; i++ ){
    //write allVars file
    allVars<<chromosome<<tab<<variantPosGlobal[i]<<endl;
  }

}

int main (int argc, char* argv[])
{

  cout.precision(13);

  ifstream fragPh(argv[1]);

  ifstream vcfFile(argv[2]);

  ofstream fragments(argv[3]);

  ofstream allVars(argv[4]);

  ofstream varData(argv[5]);

  string chromosome = argv[6];

  passHeader(vcfFile);


  //go to the right chr in the vcf files
  chargeNextChrVar(chromosome,vcfFile);

  //charge the molecules for the chr
  chargeNextChrFrag(chromosome,fragPh);

  cout <<"No of variants before coverage filter:" << variantPosGlobal.size()<<endl;

  computeMolCovPerVar();
  removeOvercovVariants();


  cout <<"No of variants after coverage filter:" << variantPosGlobal.size()<<endl;

  cout<<"Mean molecular coverage per variant:"<<meanMolCoverage<<endl;
  cout<<"Mean molecular with read coverage per variant:"<<meanMolWithReadCoverage<<endl;

  cout<<"No of variants:"<<variantPosGlobal.size()<<" "<<variantMolCoverage.size()<<" "<<
      variantMolWithReadCoverage.size()<<" "<<variantHapCoverage.size()<<endl;



  for(int i=0;i<variantPosGlobal.size();i++){

    varData<<variantPosGlobal[i]<<tab<<variantQualGlobal[i]<<tab<<
          variantMolCoverage[i]<< tab<< variantMolWithReadCoverage[i]<<tab<<
          variantHapCoverage[i].first<<tab<<variantHapCoverage[i].second<<endl;
  }
  varData.close();

  printVarsFile(chromosome, allVars);



  int noOfMolecules = allMolecules.size();

  for(int i=0; i<noOfMolecules; i++){

      printMoleculeInfo(allMolecules[i], fragments);
  }


  fragments.close();
  allVars.close();

  return 0;

}
