// -*- C++ -*-
// LHAPDFv5/v6 compatibility example

#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <string.h>

using namespace std;

LHAPDF::PDF* *pdfs;

extern "C" {

  void setlha6init_(int &maxsets) {
    pdfs = new LHAPDF::PDF* [maxsets];
  }
  
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
      if(id==0) id=21;
      fx[j] = pdfs[iset]->xfxQ2(id, x, q2);
    }
  }

  bool generic_has_id_(int &iset, int &id) {
    return pdfs[iset]->hasFlavor(id);
  }
  
  void xf_pdgid_(int &iset, int &id, double &x, double &q2, double &xf) {
      xf = pdfs[iset]->xfxQ2(id, x, q2);
  }

  void setlha6del_(int &iset) {
    delete pdfs[iset];
  }

  void alphasfrompdf0_(int &iset,double &q, double &asq) {
    asq = pdfs[iset]->alphasQ(q);
  }

  void lhapdfname_(int &ndns, char* fstring, int &iset){
    string cname=LHAPDF::lookupPDF(ndns).first;
    strcpy(fstring,cname.c_str());
    iset=LHAPDF::lookupPDF(ndns).second;
    //cout<<cname<<"|"<<endl;
    //cout<<fstring<<"|"<<endl;
  }
}
