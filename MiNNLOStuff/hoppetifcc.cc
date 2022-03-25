// -*- C++ -*-
// LHAPDFv5/v6 compatibility example

#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <string.h>

using namespace std;

LHAPDF::PDF* *pdfs;

extern "C" {


  void lhapdfname_(int &ndns, char* fstring, int &iset){
    string cname=LHAPDF::lookupPDF(ndns).first;
    strcpy(fstring,cname.c_str());
    iset=LHAPDF::lookupPDF(ndns).second;
  }

  void setlha6set_(int &iset, int &ndns, int &order, double &mz, double &asmz,double &q2min) {
    pair<string,int> set_mem = LHAPDF::lookupPDF(ndns);
    pdfs[iset] = LHAPDF::mkPDF(set_mem.first, set_mem.second);
    order = pdfs[iset]->orderQCD();
    asmz = pdfs[iset]->alphasQ(mz);
    q2min = pdfs[iset]->q2Min();
  }
}
