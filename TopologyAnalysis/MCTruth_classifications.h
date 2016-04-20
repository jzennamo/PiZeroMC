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
    int fpdg = 0;
    int nparticles = 0;
    int fProtons = 0;
      
    // X in PiZero + Mu + X neutrino interaction topolgy
    int all = 0;
    int nothing = 0;
    int nothing_n = 0;
    int single_proton = 0;
    int single_proton_n = 0;
    int single_chPi = 0;
    int single_chPi_n = 0;
    int single_pi0 = 0;
    int single_pi0_n = 0;
    int proton_pi0 = 0;
    int proton_pi0_n = 0;
    int proton_chPi = 0;
    int proton_chPi_n = 0;
    int Pi0_more2pro = 0;
    int Pi0_more2pro_n = 0;
    int proton_more2 = 0;
    int proton_more2_n = 0;
    int chPi_more2 = 0;
    int chPi_more2_n = 0;
    int chPi_more2_pros = 0;
    int chPi_more2_pros_n = 0;
    int chPi_more2pro = 0;
    int chPi_more2pro_n = 0;
    int more_chPi_pro_pi0 = 0;
    int more_chPi_pro_pi0_n = 0;
    int more_chPi_pi0 = 0;
    int more_chPi_pi0_n = 0;
    int least1Kaon = 0;
    int morePi0_pros = 0;
    int morePi0_pros_n = 0;
    int nuclear = 0;
    int photons = 0;
      
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
