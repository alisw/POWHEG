//----------------------------------------------------------------------
// example code of pythia8 shower applied on lhef files
// for the powheg directphoton process
//----------------------------------------------------------------------
// There are analysis cuts applied following an ATLAS measurement
// of isolated photons (INSPIREHEP 1510441).
// The output should reproduce fig.2 from powheg direct photon
// simulations from INSPIREHEP 1623205.
//----------------------------------------------------------------------
// author: H. Poppenborg (hendrik.popppenborg@wwu.de)
//----------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <map>
#include <cstring>
#include <vector>
#include <TStyle.h>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"

#ifdef _CPPPWHGHOOKS
#include "Pythia8Plugins/PowhegHooks.h"
#define _CPPHOOKS PowhegHooks
#else
#include "QEDQCDPowhegHooks.h"
#define _CPPHOOKS QEDQCDPowhegHooks
#endif

using namespace Pythia8;

int main(int, char **);
int main(int argc, char **argv) {

  double CorrectPhiDelta(double a, double b);
  void Fill_For_Each_Weight(vector<TH1D> &vec_h, double val, vector<double> &vec_weights);
  
  int nFiles;
  string fileName;
  Pythia p;
  _CPPHOOKS *powhegHooks = 0; // POWHEG UserHooks
  bool loadhooks;
  p.readFile("shower.conf");

  //---read commandline args----------------------------------------
  if (argc < 3) {
    cout << endl << "Usage: " << argv[0]
         << "outputfile.root eventfile1.lhe eventfile2.lhe..." << endl;
    exit(EXIT_FAILURE);
  } 
  const char *rootFileName = argv[1]; // output file
  nFiles = argc - 2; // number of event files to process

  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);

  // fill histograms with measured data points to compare with
  const int nPtBins = 14;
  const double etaAbsMin[4] = {0.0, 0.6, 1.56, 1.81};
  const double etaAbsMax[4] = {0.6, 1.37, 1.81, 2.37};
  const double ptBins[nPtBins+1] = {125.,150.,175.,200.,250.,300.,350.,400.,470.,550.,650.,750.,900.,1100.,1500.};
  
  const double ptData_rap1[14] = {2.27,1.02,0.521,0.225,0.0827,0.0357,0.0173,0.00804,0.00351,0.00133,0.000557,0.000227,0.0000925,0.0000085};
  const double ptData_rap2[14] = {2.93,1.32,0.665,0.285,0.1050,0.0452,0.0211,0.00959,0.00401,0.00166,0.000623,0.000248,0.0000569,0.0000072};
  const double ptData_rap3[13] = {0.944,0.419,0.211,0.0916,0.0334,0.0138,0.00600,0.00293,0.00110,0.000340,0.000109,0.0000297,0.0000057};
  const double ptData_rap4[12] = {1.87,0.816,0.404,0.168,0.0553,0.0217,0.00899,0.00326,0.00115,0.000317,0.000117,0.0000219};
  const double ptError_rap1[14] ={0.08,0.03,0.014,0.0058,0.0022,0.0012,0.00060,0.00031,0.00017,0.000081,0.000050,0.000026,0.000014,0.0000029};
  const double ptError_rap2[14] ={0.13,0.047,0.024,0.0093,0.0034,0.0016,0.00081,0.00040,0.00020,0.00011,0.000055,0.000028,0.000011,0.0000049};
  const double ptError_rap3[13] ={0.076,0.035,0.019,0.0081,0.0032,0.0014,0.00066,0.00035,0.00016,0.000058,0.000027,0.000013,0.0000035};
  const double ptError_rap4[12] ={0.066,0.030,0.015,0.0067,0.0025,0.0011,0.00056,0.00024,0.00011,0.000044,0.000033,0.0000078};

  vector<TH1D> vec_data;
  vector< vector<TH1D> > vec_sim; // two vector dimensions for different rap. bins & for different weights (e.g. for scale/pdf variation)

  for(int i = 0; i <= 3; i++) // book histograms for each rapidity bin
    vec_data.push_back(TH1D(Form("atlas_data_rap%d",i), Form("isolated photons in %f < |#eta| < %f;p_{T} (GeV);#frac{d#sigma}{dp_{T}}(pb)", etaAbsMin[i], etaAbsMax[i]), nPtBins, ptBins));

  for(int i = 1; i <= 14; ++i){ // fill data points
    vec_data.at(0).SetBinContent(i, ptData_rap1[i-1]);
    vec_data.at(0).SetBinError(i, ptError_rap1[i-1]);
    vec_data.at(1).SetBinContent(i, ptData_rap2[i-1]);
    vec_data.at(1).SetBinError(i, ptError_rap2[i-1]);
  }
  for(int i = 1; i <= 13; ++i){
    vec_data.at(2).SetBinContent(i, ptData_rap3[i-1]);
    vec_data.at(2).SetBinError(i, ptError_rap3[i-1]);
  }
  for(int i = 1; i <= 12; ++i){
    vec_data.at(3).SetBinContent(i, ptData_rap4[i-1]);
    vec_data.at(3).SetBinError(i, ptError_rap4[i-1]);
  }


  // prepare bookkeeping of weights
  //----------------------------------------------------------------------
  const string sudaWeightID = "sudakovwgt";
  bool isSudaWeight = false; // was photon radiation enhanced?
  double sudaWeight = 1.;    // reweighting factor associated with radiation enhancement
  vector<double> vec_weights;   // shall later contain: sudaWeight * primary event weight (using vector to store multiple weights, e.g for scale/pdf variation)
  vector<string> vec_weightsID;// vector storing descriptive id of weights

  // pythia settings required for usage with powheg
  //----------------------------------------------------------------------
  p.readString("Beams:frameType = 4");
  p.readString("Next:numberCount = 50000");

  // read in from conf file
  int vetoMode    = p.settings.mode("POWHEG:veto");
  int MPIvetoMode = p.settings.mode("POWHEG:MPIveto");
  loadhooks = (vetoMode > 0 || MPIvetoMode > 0);

  if (loadhooks) { // if NOT use SCALUP as starting scale
    if (vetoMode > 0) { // use kinematical limit as starting scale and veto
      p.readString("SpaceShower:pTmaxMatch = 2");
      p.readString("TimeShower:pTmaxMatch = 2");
    }
    if (MPIvetoMode > 0) {
      p.readString("MultipartonInteractions:pTmaxMatch = 2");
    }
    // activate POWHEG compliance
    powhegHooks = new _CPPHOOKS();
    p.setUserHooksPtr((UserHooks *) powhegHooks);
  }

  
  // variables to keep track of
  //----------------------------------------------------------------------
  int nEvents = 0,
    iPhoton = -1;     // index of photon in pythia event

  double ptMax = 0., // pT of hardest photon
    ptTemp = 0.,
    etaAbsPhoton = 999.;

  const double isoConeRadius = 0.4;
  double isoCone_et;
  double isoCone_dR;

  
  // loop over lhef files showering each event
  //----------------------------------------------------------------------
  for (int iFile = 0; iFile < nFiles; iFile++) {

    fileName = argv[2 + iFile];
    cout << "Showering events in " << fileName << endl;

    // tell Pythia to use several lhe files while initializing once
    if (iFile == 1) p.readString("Beams:newLHEFsameInit = on");
    p.readString("Beams:LHEF = " + fileName);
    p.init();

    // skip pythia errors and break, when showering has reached the end of the LHE file
    //----------------------------------------------------------------------
    while (true) {
      if (!p.next()) {
        if (p.info.atEndOfFile()) break;
        continue;
      }

      nEvents++;

      // at very first event read in weight IDs and book histograms for each weight
      //----------------------------------------------------------------------
      if (nEvents == 1 && iFile == 0) {

	// check if the sudakov weight from enhanced radiation is present
	for (map<string,double>::iterator it = p.info.weights_detailed->begin();
             it != p.info.weights_detailed->end(); ++it) {
          if (it->first == sudaWeightID.c_str()) {
            isSudaWeight = true;
            printf("Sudakov reweighting of hard process is taken into account.\n");
            continue;
          }
        }

	// // if more weights at the same time are used,
	// // e.g. for scale or pdf variation, you can  access them like this
        // for (map<string,double>::iterator it = p.info.weights_detailed->begin();
        //      it != p.info.weights_detailed->end(); ++it) {
	//   if (it->first.find("scales") != std::string::npos){
	//     vec_weightsID.push_back(it->first);
	//   }
	// }

	// insert central value always at first position for convenience 
	for (map<string,double>::iterator it = p.info.weights_detailed->begin();
	     it != p.info.weights_detailed->end(); ++it) {
	  if (it->first == "central"){ // NB: these strings follow 'lhrwgt_id' in powheg-input.save
	    vec_weightsID.insert(vec_weightsID.begin(), it->first);
	  }
	}

	printf("Number of weights = %lu\n", vec_weightsID.size());
	for(long unsigned int i = 0; i < vec_weightsID.size(); i++)
	  printf("weight description at position %lu: %s\n", i, vec_weightsID.at(i).c_str());

	// book histograms for each weight
	for(int i = 0; i <= 3; i++){ // consider each rapidity bin
	  vector <TH1D> vec_histo_temp;
	  for(long unsigned int j = 0; j < vec_weightsID.size(); j++){
	    // TH1D *p = (TH1D*)vec_data.at(i).Clone(Form("atlas_sim_rap%d_%s",i,vec_weightsID.at(j).c_str()));
	    // vec_histo_temp.push_back( &p );
	    vec_histo_temp.push_back( *(TH1D*)vec_data.at(i).Clone(Form("atlas_sim_rap%d_%s",i,vec_weightsID.at(j).c_str())) );
	  }
	  vec_sim.push_back(vec_histo_temp);
	}	 
      } // back to the general event loop...

      
      // if Sudakov reweighting is activated, get corresponding weight for this event
      if (isSudaWeight) sudaWeight = p.info.getWeightsDetailedValue(sudaWeightID.c_str());

      // reload vector with regular weights * sudaWeight for this event 
      if(vec_weights.size() != 0) vec_weights.clear();
      for(long unsigned int i = 0; i < vec_weightsID.size(); i++){
        vec_weights.push_back(p.info.getWeightsDetailedValue(vec_weightsID.at(i)) * sudaWeight);
      }

      // The actual event analysis starts here.
      ptMax  = 0.;
      ptTemp = 0.;
      iPhoton = -1;
      etaAbsPhoton = 999.;

      // search for hardest photon in this event
      //----------------------------------------------------------------------
      for (int i = 5; i < p.event.size(); i++) {
        if (p.event[i].id() == 22 && p.event[i].isFinal() && // final photon
            p.event[i].status() < 90 &&                      // no decay photons allowed, only direct photons
            TMath::Abs(p.event[i].eta()) < etaAbsMax[3] &&   // in maximal ATLAS acceptance
            p.event[i].pT() > ptBins[0]){                    // in the pt reach of interest
	  
          // find ptMax
          ptTemp = p.event[i].pT();
          if (ptTemp > ptMax) {
            ptMax = ptTemp;
            iPhoton = i; // remember index of hardest photon
          }
        }
      }

      // skip to next event, if no photon was found
      if(iPhoton < 0) continue;
      
      etaAbsPhoton = TMath::Abs(p.event[iPhoton].eta());

      // use following line to ignore events with extreme weights that can cause ugly fluctuations
      // but make sure the cross section does not decrease significantly
      if(ptMax > p.info.getScalesAttribute("uborns")*2.5){
	nEvents--;
	continue;
      }

      // isolation cut: sum energy around photon and abandon event if threshold is reached
      //----------------------------------------------------------------------
      isoCone_et = 0.; // reset sum of energy in cone
      for (int i = 5; i < p.event.size(); i++) {
	if ( !p.event[i].isFinal() ) continue;
	if ( !p.event[i].isVisible() ) continue;
	if ( TMath::Abs(p.event[i].eta()) > etaAbsMax[3]+isoConeRadius ) continue;
	if ( i == iPhoton ) continue;

	// distance between photon and particle at index i
	isoCone_dR = sqrt( pow(CorrectPhiDelta(p.event[i].phi(), p.event[iPhoton].phi()), 2)
			   + pow(p.event[i].eta() - p.event[iPhoton].eta(), 2) );

	// ATLAS isolation ignores energy close to the photon
	if( TMath::Abs(CorrectPhiDelta(p.event[i].phi() , p.event[iPhoton].phi())) > 0.0875 ||
	    TMath::Abs(p.event[i].eta() - p.event[iPhoton].eta()) > 0.0625 )
	  if(isoCone_dR < isoConeRadius) isoCone_et += p.event[i].eT();
      }

      // jump to next event if hardest photon is not isolated
      if( isoCone_et >= 4.8 + 0.0042*ptMax ) continue;
      
      // Fill histograms
      //----------------------------------------------------------------------
      for( int i = 0; i <= 3; i++)
	if( etaAbsMin[i] < etaAbsPhoton &&
	    etaAbsPhoton < etaAbsMax[i] )
	  Fill_For_Each_Weight(vec_sim.at(i), ptMax, vec_weights);
	 
    } // end of while loop; break if next file
  } // end of file loop

  // statistics on event generation
  p.stat();

  // write histograms to file ----------------------------------------
  TFile outFile(rootFileName, "RECREATE");

  // normalize simulated spectra for nEvents and pt bin width, then write
  for(int i = 0; i <= 3; i++){
    vec_data.at(i).Write();
    for(unsigned long int j = 0; j < vec_weights.size(); j++){
      vec_sim.at(i).at(j).Scale( 1./nEvents, "width");
      vec_sim.at(i).at(j).Write();
    }
  }
  
  outFile.Close();

  if (powhegHooks) delete powhegHooks;
  return 0;

}

// PYTHIA8's phi goes from -pi to pi; compute correct angle difference
//----------------------------------------------------------------------
double CorrectPhiDelta(double angle1, double angle2){
  double pi = TMath::Pi();
  double phi = TMath::Abs(angle1 - angle2);
  if(phi >= pi) return 2*pi-phi;
  else return phi;
}

//----------------------------------------------------------------------
void Fill_For_Each_Weight(vector<TH1D> &vec_h, double val, vector<double> &vec_weights){

  for(unsigned long int i = 0; i < vec_weights.size(); i++)
    vec_h.at(i).Fill(val,vec_weights.at(i));

  return;
}

//----------------------------------------------------------------------
