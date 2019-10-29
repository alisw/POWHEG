// -*- C++ -*-
//
// powhegAnalysis.h - (c) Silvia Ferrario Ravasio and Tomas Jezo
// inspired by HepMCFile.h which is a part of ThePEG
//

#include "powhegAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/HepMCHelper.h"
#include "HepMC/IO_GenEvent.h"

#include "herwig7_interface.h"

using namespace ThePEG;
extern int maxev;
bool showerSuccess;

powhegAnalysis::powhegAnalysis() 
  : _runNumber(0), _unitchoice(),
    _geneventPrecision(16) {}

// Cannot copy streams. 
// Let doinitrun() take care of their initialization.
powhegAnalysis::powhegAnalysis(const powhegAnalysis & x) 
  : AnalysisHandler(x), _runNumber(x._runNumber),
    _unitchoice(x._unitchoice), _geneventPrecision(x._geneventPrecision) {}

IBPtr powhegAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr powhegAnalysis::fullclone() const {
  return new_ptr(*this);
}

void powhegAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();	  
  hepevt_.nevhep=0;
}

void powhegAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  herwig7_end_(&(_runNumber));
  cout << "\npowhegAnalysis: analysis finished.\n";
  cout << "Number of analyzed events: "<<hepevt_.nevhep<<endl;
  cout.flush();
}

void powhegAnalysis::analyze(tEventPtr event, long, int, int) {
  showerSuccess=true;

  Energy eUnit;
  Length lUnit;
  switch (_unitchoice) {
  default: eUnit = GeV; lUnit = millimeter; break;
  case 1:  eUnit = MeV; lUnit = millimeter; break;
  case 2:  eUnit = GeV; lUnit = centimeter; break;
  case 3:  eUnit = MeV; lUnit = centimeter; break;
  }

  HepMC::GenEvent * hepmc 
    = HepMCConverter<HepMC::GenEvent>::convert(*event, false,
					       eUnit, lUnit);

  if( hepmc->is_valid() )
    {
      int jhep = 0; //particles counter
      int maxjhep=0;
		
		
      //Loop over the particles
      for ( HepMC::GenEvent::particle_const_iterator p =
	      hepmc->particles_begin(); p != hepmc->particles_end(); ++p)
	{
	  //(**p).print();
	  //Copy the properties of the particles in the hepevt 
	  //common block. (HepMC barcode starts from 10001, our 
	  //fortran array will start from 1).
	  jhep  = (*p)->barcode() - 10001; //C++ array starts from 0
	  if(jhep>maxjhep) maxjhep = jhep;
			
	  hepevt_.idhep[jhep]   = (*p)->pdg_id();
	  hepevt_.isthep[jhep]  = (*p)->status();
	  hepevt_.phep[jhep][0] = (*p)->momentum().px()	;
	  hepevt_.phep[jhep][1] = (*p)->momentum().py();
	  hepevt_.phep[jhep][2] = (*p)->momentum().pz();
	  hepevt_.phep[jhep][3] = (*p)->momentum().e();
	  hepevt_.phep[jhep][4] = (*p)->momentum().m();			

	  //parents location
	  if ( (*p)->production_vertex() ) //the particles has 1 or 2 moms..
	    {
	      int nmom =0;  
	      for ( HepMC::GenVertex::particle_iterator mother = 
		      (*p)->production_vertex()->
		      particles_begin(HepMC::parents);
		    mother != (*p)->production_vertex()->
		      particles_end(HepMC::parents); ++mother )
		{
		  hepevt_.jmohep[jhep][nmom]=(*mother)->barcode()-10000; // this is assigning the fortran index
		                                                         // So -10000 and NOT -10001
		  nmom ++;
		}
	      if (nmom == 1) 
		hepevt_.jmohep[jhep][1]=hepevt_.jmohep[jhep][0];		      			     
	    }
	  else  //incoming protons: 0 moms
	    {
	      hepevt_.jmohep[jhep][0]=0;
	      hepevt_.jmohep[jhep][1]=0;
	    }

	  // children location. We need the first and the last
	  if ( (*p)->end_vertex() ) 
	    {
	      int nchild =0;
	      //std::cout<<"children: ";
	      for ( HepMC::GenVertex::particle_iterator child 
		      =(*p)->end_vertex()->
		      particles_begin(HepMC::children);
		    child != (*p)->end_vertex()->
		      particles_end(HepMC::children);++child ) 
		{
		  if(nchild == 0) hepevt_.jdahep[jhep][0]=(*child)->barcode()-10000; //first son
		  hepevt_.jdahep[jhep][1]=(*child)->barcode()-10000; //override the last son
		  nchild ++; //number of sons
		}
	    }
	  else //final state particles
	    {
	      hepevt_.jdahep[jhep][0] = -1;
	      hepevt_.jdahep[jhep][1] = -1;
	    }
	}
      hepevt_.nhep= maxjhep +1;
    }
     	

  hepevt_.nevhep++;
  herwiganalysis_();
  if( (hepevt_.nevhep % 2000)==0) herwig7_end_(&(_runNumber));
	

  delete hepmc;
  if(hepevt_.nevhep>=maxev)
    {
      powhegAnalysis::dofinish();
      exit(0);
      return;      
    }
}

void powhegAnalysis::persistentOutput(PersistentOStream & os) const {
  os << _runNumber << _unitchoice << _geneventPrecision;
}

void powhegAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> _runNumber >> _unitchoice >> _geneventPrecision;
}

ClassDescription<powhegAnalysis> powhegAnalysis::initpowhegAnalysis;
// Definition of the static class description member.

void powhegAnalysis::Init() {

  static ClassDocumentation<powhegAnalysis> documentation
    ("This analysis handler will convert first into HepMC format then"
     "into hepevt common block and then run the analysis implemented in"
     "the file pwhg_analysis*.f (specified in the Makefile).");

  static Parameter<powhegAnalysis,int> interfaceRunNumber
    ("RunNumber",
     "The number identifying the run. The run number will be used in the filename of the"
     ".top file produced. For example, if run number is 2, `pwg-0002-POWHEG+HERWIG7-output.top`"
     "file containing the histograms will be produced. Values below 1 and above 9999 will be"
     "ignored.",
     &powhegAnalysis::_runNumber,
     0, // the default value
     1, // the minimum
     9999, // the maximum
     true, // depsafe
     false, // readonly
     Interface::limited); //limits

  static Parameter<powhegAnalysis,unsigned int> interfacePrecision
    ("Precision",
     "Choice of output precision for the HepMC GenEvent format "
     " (as number of digits).",
     &powhegAnalysis::_geneventPrecision, 16, 6, 16,
     false, false, Interface::limited);
  
  static Switch<powhegAnalysis,int> interfaceUnits
    ("Units",
     "Unit choice for energy and length",
     &powhegAnalysis::_unitchoice, 0, false, false);
  static SwitchOption interfaceUnitsGeV_mm
    (interfaceUnits,
     "GeV_mm",
     "Use GeV and mm as units.",
     0);
  static SwitchOption interfaceUnitsMeV_mm
    (interfaceUnits,
     "MeV_mm",
     "Use MeV and mm as units.",
     1);
  static SwitchOption interfaceUnitsGeV_cm
    (interfaceUnits,
     "GeV_cm",
     "Use GeV and cm as units.",
     2);
  static SwitchOption interfaceUnitsMeV_cm
    (interfaceUnits,
     "MeV_cm",
     "Use MeV and cm as units.",
     3);
}
