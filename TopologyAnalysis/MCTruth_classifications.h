/**
 * \file MCTruth_classifications.h
 *
 * \ingroup TopologyAnalysis
 * 
 * \brief Class def header for a class MCTruth_classifications
 *
 * @author jzennamo
 */

/** \addtogroup TopologyAnalysis

    @{*/

#ifndef LARLITE_MCTRUTH_CLASSIFICATIONS_H
#define LARLITE_MCTRUTH_CLASSIFICATIONS_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcpart.h"
#include "TTree.h"

namespace larlite {
  /**
     \class MCTruth_classifications
     User custom analysis class made by SHELL_USER_NAME
   */
  class MCTruth_classifications : public ana_base{
  
  public:

    /// Default constructor
    MCTruth_classifications(){ _name="MCTruth_classifications"; _fout=0;}

    /// Default destructor
    virtual ~MCTruth_classifications(){}

    /** IMPLEMENT in MCTruth_classifications.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MCTruth_classifications.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MCTruth_classifications.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    // Variables you want
    int all_events;
    int CC_events;
    int CCpi0_events; 
    
    //Variable for your tree go here!
    TTree* tree;
    int thisInt = 0;
    double thisDouble = 0.0;

  protected:
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
