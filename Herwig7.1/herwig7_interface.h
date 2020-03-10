// herwig7_interface.h - (c) Silvia Ferrario Ravasio, Tomas Jezo, 
//   Paolo Nason and Carlo Oleari

#define max_num_weights  50

extern "C" {
  void herwiganalysis_();
  void herwig7_end_(int*);
  void herwig7_init_(int & maxev, const char* name, int len);
  double powheginput_(const char* , int);

  //Declaration of fortran common blocks

  extern struct {
    int idbmup[2];
    double ebmup[2];
    int pdfgup[2], pdfsup[2], idwtup, nprup;
    double xsecup[100], xerrup[100], xmaxup[100];
    int lprup[100];
  } heprup_;

  extern struct {
    int nup, idprup;
    double xwgtup, scalup, aqedup, aqcdup;
    int idup[500], istup[500], mothup[500][2], icolup[500][2];
    double pup[500][5], vtimup[500],spinup[500];
  } hepeup_;

  extern struct {	  
    int nevhep, nhep, isthep[4000], idhep[4000], jmohep[4000][2], 
      jdahep[4000][2];
    double phep[4000][5], vhep[4000][4];
  } hepevt_;

  extern struct{
    double weight[max_num_weights];
    int numweights;
    int radtype;
  } weights_;
  
}

