#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;


int main (int argc, char* argv[])
{

  string interleavedFastq=argv[1];
  ifstream interFastq(interleavedFastq.c_str());

  string filtered_interleavedFastq=argv[2];
  ofstream filt_interFastq(filtered_interleavedFastq.c_str());

  string name, sequence, plus, quality, barcode, commentBC;
  string prevName="";
  string pair="";

  int bcPos=0;

  int totalReads=0;
  int pairedReads=0;

  while(getline(interFastq, name)){

    totalReads++;

    if(prevName.compare(name) != 0){
      pair=":1";
    }else{
      pair=":2";
      pairedReads++;
    }

    prevName=name;

    getline(interFastq, sequence);
    getline(interFastq, plus);
    getline(interFastq, quality);

    bcPos = name.find(" BX:Z:");

    if( bcPos != string::npos){

      barcode = name.substr(bcPos+5);
      commentBC = name.substr(bcPos);
      name = name.substr(0,bcPos);
      //name.append(pair);
      name.append(barcode);
      name.append(commentBC);


      filt_interFastq << name << endl;
      filt_interFastq << sequence << endl;
      filt_interFastq << plus << endl;
      filt_interFastq << quality << endl;

    }

  }

  filt_interFastq.close();

  cout<<totalReads<<":"<<pairedReads;

  return 0;
}
