/**
 * \file MCTrack_lengths.h
 *
 * \ingroup TopologyAnalysis
 * 
 * \brief Class def header for a class MCTrack_lengths
 *
 * @author laube
 */

/** \addtogroup TopologyAnalysis

    @{*/

#ifndef LARLITE_MCTRACK_LENGTHS_H
#define LARLITE_MCTRACK_LENGTHS_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
//#include "UserDev/BasicTool/GeoAlgo/GeoAABox.h"
//#include "GeoAlgo/GeoAlgo.h"
#include "TTree.h"

namespace larlite {
  /**
     \class MCTrack_lengths
     User custom analysis class made by SHELL_USER_NAME
   */
  class MCTrack_lengths : public ana_base{
  
  public:

    /// Default constructor
    MCTrack_lengths(){ _name="MCTrack_lengths"; _fout=0;}

    /// Default destructor
    virtual ~MCTrack_lengths(){}

    /** IMPLEMENT in MCTrack_lengths.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MCTrack_lengths.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MCTrack_lengths.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

      //Variables you want
      
      //Variables for the tree go here!
      
      TTree* tree;
      int thisInt           = 0;
      int no_muons          = 0;
      int no_tracks         = 0;
      int no_prim_tracks    = 0;
      
      TTree* muons;
      double muon_length    = 0;
      
      TTree* other;
      double other_length   = 0;
      double other_pdg      = 0;
      std::string other_process;
      double d_vtx          = 0;
      
      int CCpi0_events      = 0;
      int all               = 0;
      int nTopoCut          = 0;
      
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
