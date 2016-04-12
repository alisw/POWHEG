// -*- C++ -*-
// LHAPDFv5/v6 compatibility example

#include "LHAPDF/LHAPDF.h"
#include <iostream>

using namespace std;


LHAPDF::PDF* pdfs[10];

extern "C" {
  void setlha6set_(int &iset, int &ndns, int &order, double &mz, double &asmz,double &q2min) {
    pair<string,int> set_mem = LHAPDF::lookupPDF(ndns);
    pdfs[iset] = LHAPDF::mkPDF(set_mem.first, set_mem.second);
    order = pdfs[iset]->orderQCD();
    asmz = pdfs[iset]->alphasQ(mz);
    q2min = pdfs[iset]->q2Min();
  }
  

  void xfxq2_(int &iset, double &x, double &q2, double fx[13]) {
    int j,id;
    for(j=0;j<13;j++) {
      id = j-6;
      if(id==0) id==21;
      fx[j] = pdfs[iset]->xfxQ2(id, x, q2);
    }
  }
  void xfphoton_(int &iset, double &x, double &q2, double photon) {
      photon = pdfs[iset]->xfxQ2(22, x, q2);
  }
}
