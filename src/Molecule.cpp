#include "Molecule.h"
#include <iostream>
using namespace std;

Molecule::Molecule(string chr, int start, int end, string bc, int noReads) :
chromosome(chr), start(start), end(end), barcode(bc), noReads(noReads)
{}

string Molecule::getChr() const{
  return this->chromosome;
}

int Molecule::getStart() const{
  return this->start;
}

int Molecule::getEnd() const{
  return this->end;
}

string Molecule::getBC() const{
  return this->barcode;
}

int Molecule::getNoReads() const{
  return this->noReads;
}

bool Molecule::spanVar(int posVar){

  bool span=false;

  if ( (this->start <= posVar) && (this->end >= posVar) ){
    span=true;
  }

  return span;
}
