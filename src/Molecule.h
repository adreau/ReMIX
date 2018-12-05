#ifndef MOLECULE_H
#define MOLECULE_H

#include <iostream>
#include <string>
using namespace std;

class Molecule
{
public:

  Molecule(string chr, int start, int end, string bc, int noReads);

  string getChr() const;
  int getStart() const;
  int getEnd() const;
  string getBC() const;
  int getNoReads() const;

  bool spanVar(int posVar);


private:

  string chromosome;
  int start, end;
  string barcode;
  int noReads;


};

#endif
