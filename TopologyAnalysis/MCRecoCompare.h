/**
 * \file MCRecoCompare.h
 *
 * \ingroup TopologyAnalysis
 * 
 * \brief Class def header for a class MCRecoCompare
 *
 * @author laube
 */

/** \addtogroup TopologyAnalysis

    @{*/

#ifndef LARLITE_MCRECOCOMPARE_H
#define LARLITE_MCRECOCOMPARE_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/pfpart.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
//#include "GeoAlgo/GeoAABox.h"
#include "LArUtil/Geometry.h"
#include "TTree.h"
#include <TH1.h>

namespace larlite {
  /**
     \class MCRecoCompare
     User custom analysis class made by SHELL_USER_NAME
   */
  class MCRecoCompare : public ana_base{
  
  public:

    /// Default constructor
    MCRecoCompare(){ _name="MCRecoCompare"; _fout=0;}

    /// Default destructor
    virtual ~MCRecoCompare(){}

    /** IMPLEMENT in MCRecoCompare.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in MCRecoCompare.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in MCRecoCompare.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();
      
      //Variables you want
      
      //Variables for the tree go here!
      
      TTree* tree;
      int thisInt           = 0;
      double lengthRatio    = 0;
      double vertexdiff     = 0;

      
      TTree* mctree;
      int mcevent           = -1;
      double muon_length    = -1;
      
      TTree* recotree;
      int pfpevent          = -1;
      double d_vtx_trk      = 0; //track distances from vertex
      double d_vtx_show      = 0; //shower distances from vertex
      
      TTree* fake;
      
      int all_events        = 0;
      int CC_events      = 0;
      int numuCCpi0_events               = 0;
      int nTopoCut          = 0;
      
      int p_CCpi0_events    = 0;
      int p_all             = 0;
      int p_nTopoCut        = 0;
      
      int nSuccess          = 0;
      int nFake             = 0;
      int fakeNotCC         = 0; // number of fakes that are not CC in MC
      int fakeNotTopo       = 0; // number of fakes that do not have the required topology (number of showers, tracks)
      int nMoreThanOne      = 0;

      int evt_id            = 0; // only for debug
      int nCCpi0Pand        = 0;
      
      int fakeNotTopoEvent = 0;
      
      TH1F *mc_track_lengths    = new TH1F("mc_track_lengths", "Lengths of tracks from vertex in selected MC events", 100, -10, 1000);
      TH1F *mc_muon_lengths     = new TH1F("mc_muon_lengths", "Lengths of muon tracks from interaction", 100, -10, 1000);
      TH1F *mc_ntracks          = new TH1F("mc_ntracks", "Number of tracks from interaction vertex", 10, 0, 10);
      
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
      
      bool StartInTPC(larlite::mcpart mctrk){
          if(mctrk.Trajectory().at(0).Position().X() < 0 || mctrk.Trajectory().at(0).Position().X() > 2*(larutil::Geometry::GetME()->DetHalfWidth()) ||
             mctrk.Trajectory().at(0).Position().Y() < -(larutil::Geometry::GetME()->DetHalfHeight()) || mctrk.Trajectory().at(0).Position().Y() > larutil::Geometry::GetME()->DetHalfHeight() ||
             mctrk.Trajectory().at(0).Position().Z() < 0 || mctrk.Trajectory().at(0).Position().Z() > larutil::Geometry::GetME()->DetLength()) {
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
