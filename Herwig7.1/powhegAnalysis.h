// -*- C++ -*-
//
// powhegAnalysis.h - (c) Silvia Ferrario Ravasio and Tomas Jezo
// inspired by HepMCFile.h which is a part of ThePEG
//
#ifndef powhegAnalysis_H
#define powhegAnalysis_H
//
// This is the declaration of the powhegAnalysis class.
//
#include <iostream>
#include <fstream>
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "HepMC/IO_BaseClass.h"

namespace ThePEG {

/** \ingroup Analysis
 * The powhegAnalysis class outputs ThePEG events in HepMC format.
 *
 * @see \ref powhegAnalysisInterfaces "The interfaces"
 * defined for powhegAnalysis.
 */
class powhegAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  powhegAnalysis();

  /**
   * The copy constructor.
   */
  powhegAnalysis(const powhegAnalysis &);
  //@}

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);
  //@}

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
   virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<powhegAnalysis> initpowhegAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  powhegAnalysis & operator=(const powhegAnalysis &);

private:

  /**
   *  The run number
   *  The number identifying the run. The run number will be used in the filename of the
   *  .top file produced. For example, if run number is 2, `pwg-0002-POWHEG+HERWIG7-output.top`
   *  file containing the histograms will be produced. Values below 1 and above 9999 will be
   *  ignored.
   */
  int _runNumber;

  /**
   * Selector for the choice of units
   */
  int _unitchoice;

  /**
   * Choice of output precision in GenEvent format
   */
  unsigned int _geneventPrecision;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of powhegAnalysis. */
template <>
struct BaseClassTrait<powhegAnalysis,1> {
  /** Typedef of the first base class of powhegAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the powhegAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<powhegAnalysis>
  : public ClassTraitsBase<powhegAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::powhegAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the powhegAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "powhegHerwig.so"; }
};

/** @endcond */

}

#endif /* THEPEG_powhegAnalysis_H */
