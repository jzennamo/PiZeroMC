/**
 * \file numuCCfilterstudy.h
 *
 * \ingroup TopologyAnalysis
 * 
 * \brief Class def header for a class numuCCfilterstudy
 *
 * @author laube
 */

/** \addtogroup TopologyAnalysis

    @{*/

#ifndef LARLITE_NUMUCCFILTERSTUDY_H
#define LARLITE_NUMUCCFILTERSTUDY_H

#include "Analysis/ana_base.h"
#include "DataFormat/mcpart.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"
#include "LArUtil/Geometry.h"
#include <TH1.h>

namespace larlite {
  /**
     \class numuCCfilterstudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class numuCCfilterstudy : public ana_base{
  
  public:

    /// Default constructor
    numuCCfilterstudy(){ _name="numuCCfilterstudy"; _fout=0;}

    /// Default destructor
    virtual ~numuCCfilterstudy(){}

    /** IMPLEMENT in numuCCfilterstudy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in numuCCfilterstudy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in numuCCfilterstudy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();
      
      
      int all_events        = 0;
      int numuCCpi0_events  = 0;
      int nTopoCut          = 0;
      int twoPi0            = 0;
      int noPi0             = 0;
      
      double muon_length    = -1;
      TVector3 muontrackstart;
      TVector3 reco_muontrackstart;
      TVector3 reco_muontrackend;
      
      TH1F *mc_track_lengths    = new TH1F("mc_track_lengths", "Lengths of tracks from vertex in selected MC events", 100, -10, 1000);
      TH1F *mc_muon_lengths     = new TH1F("mc_muon_lengths", "Lengths of muon tracks from interaction", 100, 0, 1000);
      TH1F *reco_muon_lengths   = new TH1F("reco_muon_lengths", "Length of reco muon candidate selected by numuCCInc filter", 100, 0, 0);
      TH1F *diff_muon_length    = new TH1F("diff_muon_length", "Difference in muon track length MC - Reco", 100, -700, 300);
      TH1F *mc_ntracks          = new TH1F("mc_ntracks", "Number of tracks from interaction vertex", 15, 0, 15);
      TH1F *vtx_distance        = new TH1F("vtx_distance", "Distance between true vertex and NuMuCCInc selected", 100, 0, 350);
      TH1F *reco_ntracks        = new TH1F("reco_ntracks", "Number of tracks reconstructed as neutrino daughters" , 15, 0, 15);
      TH1F *mc_nshowers         = new TH1F("mc_nshowers", "Number of true photon showers", 6, 0, 6);
      TH1F *reco_nshowers       = new TH1F("reco_nshowers", "Number of reconstructed showers from vertex", 6, 0, 6);
      TH1F *diff_ntracks        = new TH1F("diff_ntracks", "Difference in number of tracks MC-Reco", 17, -5, 12 );
      TH1F *diff_nshowers       = new TH1F("diff_nshowers", "Difference in number of showers MC-Reco", 6, -3, 3);
      TH1F *diff_muon_length_cut= new TH1F("diff_muon_length_cut", "Difference in muon track length for muons > 75cm, MC-Reco", 100, -700, 300);
      TH1F *vtx_distance_cut    = new TH1F("vtx_distance_cut", "Distance between true vertex and NuMuCCInc selected for true muons < 75cm", 100, 0, 350);
      TH1F *diff_nclosetracks   = new TH1F("diff_nclosetracks", "Difference in number of close to vertex tracks MC-Reco", 17, -5, 12);
      TH1F *TrackVertexDist     = new TH1F("TrackVertexDist", "Distance of Start of Nu Daughter Track to Nu Vertex", 100, 0,0);
      TH1F *dist_muonStart      = new TH1F("dist_muonStart", "Distance between MC and Reco Muon Start", 100, 0, 0);
      TH1F *dist_muonStart_cut  = new TH1F("dist_muonStart_cut", "Distance between MC and Reco Muon Start tracks >75cm", 100, 0, 0);
      TH1F *total_ntrk_nshw     = new TH1F("total_ntrk_nshw", "Difference in (Number of tracks + number of showers) MC - Reco", 30,-15,15);
      
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
      
      bool _verbose;
      int _n_tracks;
      int _n_showers;
      
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
//
//
//**************************************************************************

/** @} */ // end of doxygen group 
