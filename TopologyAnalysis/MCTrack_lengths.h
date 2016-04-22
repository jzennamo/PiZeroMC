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
//#include "GeoAlgo/GeoAABox.h"
#include "LArUtil/Geometry.h"
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
      double d_vtx_trk      = 0; //track distances from vertex
      double d_vtx_show      = 0; //shower distances from vertex
      
      int CCpi0_events      = 0;
      int all               = 0;
      int nTopoCut          = 0;
      
  private:
      
      bool StartInTPC(larlite::mctrack mctrk){
          if(mctrk.Start().X() < 0 || mctrk.Start().X() > 2*(larutil::Geometry::GetME()->DetHalfWidth()) ||
             mctrk.Start().Y() < -(larutil::Geometry::GetME()->DetHalfHeight()) || mctrk.Start().Y() > larutil::Geometry::GetME()->DetHalfHeight() ||
             mctrk.Start().Z() < 0 || mctrk.Start().Z() > larutil::Geometry::GetME()->DetLength()) {
              return false;
          }
          else return true;
      }
      
      bool StartInTPC(larlite::mcshower mctrk){
          if(mctrk.Start().X() < 0 || mctrk.Start().X() > 2*(larutil::Geometry::GetME()->DetHalfWidth()) ||
             mctrk.Start().Y() < -(larutil::Geometry::GetME()->DetHalfHeight()) || mctrk.Start().Y() > larutil::Geometry::GetME()->DetHalfHeight() ||
             mctrk.Start().Z() < 0 || mctrk.Start().Z() > larutil::Geometry::GetME()->DetLength()) {
              return false;
          }
          else return true;
      }
      
      double length(larlite::mctrack mctrk) {
          double len = 0;
          bool first = true;
          TVector3 disp;
          for(auto step : mctrk) {
              if(first) {
                  first = false;
                  disp = {step.X(), step.Y(), step.Z()};
              }
              else {
                  TVector3 pos(step.X(), step.Y(), step.Z());
                  disp -= pos;
                  len += disp.Mag();
                  disp = pos;
              }
          }
          return len;
      }
      
      
             
             
      
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
