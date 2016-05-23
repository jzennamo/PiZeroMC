#ifndef LARLITE_NUMUCCFILTERSTUDY_CXX
#define LARLITE_NUMUCCFILTERSTUDY_CXX

#include "numuCCfilterstudy.h"
#include "DataFormat/event_ass.h"
#include "DataFormat/track.h"
#include "DataFormat/vertex.h"
#include "DataFormat/pfpart.h"



namespace larlite {

  bool numuCCfilterstudy::initialize() {
      
      _verbose = true;
      _n_tracks = 0;
      _n_showers = 0;

  // General
      
      all_events        = 0;
      numuCCpi0_events  = 0;
      nTopoCut          = 0;
      twoPi0            = 0;
      noPi0             = 0;

    
    return true;
  }
  
  // ===========================================================================================//
    
    // Manual:
    // http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
    //
    
    
  bool numuCCfilterstudy::analyze(storage_manager* storage) {
      
      
      // print event information
      std::cout << "-----------------------------------------" << std::endl;
      std::cout << "Run   " << storage->run_id() << std::endl;
      std::cout << "Event " << storage->event_id() << std::endl;
      all_events++;
      
      
      
      //check whether there is a true PiZero in Event
      
      // Neutrino Vertex
      TVector3 vtx_nu = {-9999,-9999,-9999};
      // Pion Vertex
      TVector3 vtx_pi = {-9999,-9999,-9999};
      
      int nPi0          = 0;
      int nMu           = 0;
      bool numuCCpi0    = false;
      bool topoCut      = false;
      
      auto mc = storage->get_data<event_mctruth>("generator");      // these are MCParticles from generator level
      
      // get the neutrino
      auto& nu = mc->at(0).GetNeutrino();
      std::cout << "*Generator neutrino " << nu.Nu().PdgCode() << std::endl;
      vtx_nu.SetXYZ(nu.Nu().Trajectory().at(0).Position().X(),nu.Nu().Trajectory().at(0).Position().Y(),nu.Nu().Trajectory().at(0).Position().Z());
      std::cout << "*neutrino vertex \t " << vtx_nu.X() << ", " << vtx_nu.Y() << ", " << vtx_nu.Z() << std::endl;
      for(auto d: nu.Nu().Daughters()) {
          std::cout << "MC Neutrino daughters: " << d <<std::endl;
      }
      
      //get particles
      auto& all_particles = mc->at(0).GetParticles();
      std::cout << "\n" << std::endl;
      for(auto par : all_particles) {
          //get PiZero
          if(par.PdgCode() == 111 && par.StatusCode() == 1 && StartInTPC(par)) {
          nPi0++;
          std::cout << "Generator PiZero " << par.PdgCode()  << " \t Status Code: " << par.StatusCode() << std::endl;
          vtx_pi.SetXYZ(par.Trajectory().at(0).Position().X(),par.Trajectory().at(0).Position().Y(),par.Trajectory().at(0).Position().Z());
          std::cout << "PiZero Vertex \t\t" << vtx_pi.X() << ", " << vtx_pi.Y() << ", " << vtx_pi.Z() << std::endl;
          }
      }
      
      // check Muon
      if(nPi0 == 1) {
          for(auto const& par : all_particles) {
              if(abs(par.PdgCode()) == 13 && StartInTPC(par) && par.Mother() == 0)  {
                  TVector3 muonvtx;
                  muonvtx.SetXYZ(par.Trajectory().at(0).Position().X(),par.Trajectory().at(0).Position().Y(),par.Trajectory().at(0).Position().Z());
                  if((muonvtx - vtx_nu).Mag() < 0.01) nMu++;
              }
          }
      }
      
      if(nPi0 == 1 && nMu == 1) {
          numuCCpi0_events++;
          numuCCpi0 = true;
      }
      
      if(nPi0 >= 2) twoPi0++;
      if(nPi0 == 0) noPi0++;
      
      
      // Get mctracks and mcshowers
      auto mctrk_v = storage->get_data<larlite::event_mctrack>("mcreco");        // this is a vector of MCtracks
      //std::cout << "mctrk_v length " << mctrk_v->size()<< std::endl;
      auto mcshow_v = storage->get_data<event_mcshower>("mcreco");      // this is a vector of MCshowers
      //std::cout << "mcshow_v length " << mcshow_v->size() << std::endl;
      
      // assess what event looks like in detector: use MCTracks and MCShowers
      std::cout << "\n MCTrack and MCShower profile" << std::endl;
      
      int n_s_photon  = 0;        // number of photon showers from Pi0 // NEED TO BE EXACTLY TWO
      int n_m         = 0;        // number of muons from vertex // NEEDS TO BE OBSERVABLE
      int n_t         = 0;        // number of tracks from vertex // COUNT THESE AND SPLIT UP COUNT

      
    if(numuCCpi0) {
    
        //count tracks
      for(auto mctrk : *mctrk_v) {
            // muon mctracks
          if(abs(mctrk.PdgCode()) == 13 && mctrk.Process() == "primary") {
            muontrackstart.SetXYZ(mctrk.Start().X(), mctrk.Start().Y(), mctrk.Start().Z());
              if(length(mctrk) > 0.5 && (muontrackstart - vtx_nu).Mag() < 0.01) {
                  n_m++;
                  muon_length = length(mctrk);
                  std::cout << "*MCTrack muon distance to vtx " << (muontrackstart - vtx_nu).Mag() << "\t len: " << muon_length << std::endl;
              }
          }
          // other mctracks
          if(abs(mctrk.PdgCode()) != 13 && mctrk.PdgCode() != 2112 && mctrk.PdgCode() != 130 && mctrk.PdgCode() != 310 && mctrk.PdgCode() != 111
             && mctrk.PdgCode() != 311 && mctrk.PdgCode() < 10000 && mctrk.Process() == "primary"  && StartInTPC(mctrk) &&length(mctrk) > 0.5) {
              TVector3 trackstart;
              trackstart.SetXYZ(mctrk.Start().X(), mctrk.Start().Y(), mctrk.Start().Z());
              if((trackstart - vtx_nu).Mag() < 0.01 && length(mctrk) > 0.5) {
                  n_t++;
                  std::cout << "*MCTrack " << mctrk.PdgCode() << " distance to vtx " << (trackstart - vtx_nu).Mag() << "\t len: " << length(mctrk) << std::endl;
              }
          }
      }
          
      //count showers
      for(auto mcshw : *mcshow_v) {
          if(mcshw.PdgCode() == 22 && StartInTPC(mcshw)) { // && mcshw.MotherPdgCode() == 111) {                  //count photon showers from pizero
                TVector3 svtx;
                svtx.SetXYZ(mcshw.Start().X(), mcshw.Start().Y(), mcshw.Start().Z());
                if((svtx - vtx_nu).Mag() < 0.01) {
                    n_s_photon++;
                    std::cout << "**MCShower Photon from Pi0 distance to vtx : " << (svtx - vtx_nu).Mag() << "\t Mother Pdg: " << mcshw.MotherPdgCode() << std::endl;
                }
          }
        }
    
      if(n_s_photon >= 2 && n_m == 1 && n_t >=0) {                            // selection Criteria (remove other shower cut)
          topoCut = true;
          nTopoCut++;
      }
          
      if(topoCut) {
          mc_ntracks->Fill(n_t + n_m);
          mc_nshowers->Fill(n_s_photon);
          mc_muon_lengths->Fill(muon_length);
          mc_track_lengths->Fill(muon_length);
          for(auto mctrk: *mctrk_v) {
              if(abs(mctrk.PdgCode()) != 13 && mctrk.PdgCode() != 2112 && mctrk.PdgCode() != 130 && mctrk.PdgCode() != 310 && mctrk.PdgCode() != 111
                  && mctrk.PdgCode() != 311 && mctrk.PdgCode() < 10000 && mctrk.Process() == "primary"  && StartInTPC(mctrk) &&length(mctrk) > 0.5) {
                  TVector3 trackstart;
                  trackstart.SetXYZ(mctrk.Start().X(), mctrk.Start().Y(), mctrk.Start().Z());
                  if((trackstart - vtx_nu).Mag() < 0.01 && length(mctrk) > 0.5) {                                                // count other tracks from vertex
                      mc_track_lengths->Fill(length(mctrk));
                  }
              }
        
          }
      }
        
    }
    // ------------ reco -------------------------------------------------------------------------------------------------------------
      
      if(!topoCut) std::cout << "\n MC topoCut false" << std::endl;
      if(topoCut) {
          std::cout << "\n MC topoCut true \n " << std::endl;
      
          
          // keep track of number of showers and tracks found
          int n_showers     = 0;
          int n_tracks      = 0;
          int n_close_track = 0;
          // keep track of the reconstructed tracks produced in the neutrino
          // interaction. This vector contains all reco tracks with parent
          // the reconstructed neutrino
          std::vector<larlite::track> nu_trk_v;

          // get a handle to the association
          auto ev_ass = storage->get_data<larlite::event_ass>("NuMuCCInclusive");
          
          // get the association keys
          auto const& ass_keys = ev_ass->association_keys();
          
          larlite::AssSet_t ass_trk_vtx_v;
          larlite::event_track *ev_trk = nullptr;
          ass_trk_vtx_v = storage->find_one_ass( ass_keys[0].second, ev_trk, ev_ass->name() );
          
          larlite::AssSet_t ass_vtx_trk_v;
          larlite::event_vertex *ev_vtx = nullptr;
          ass_vtx_trk_v = storage->find_one_ass( ass_keys[0].first, ev_vtx, ev_ass->name() );
          
          // are there tracks? are there vertices?
          if (!ev_trk or (ev_trk->size() == 0)){
              std::cout << "No track! exit" << std::endl;
              return false;
          }
          if (!ev_vtx or (ev_vtx->size() == 0)){
              std::cout << "No vertex! exit" << std::endl;
              return false;
          }
          
          // grab PFParticles associated with these tracks
          larlite::AssSet_t ass_trk_pfpart_v;
          larlite::event_pfpart *ev_pfpart = nullptr;
          bool pfpart = true;  // have the PFParticles been found?
          ass_trk_pfpart_v = storage->find_one_ass( ev_trk->id(), ev_pfpart, ev_trk->name() );
          
          if (!ev_pfpart or (ev_pfpart->size() == 0)){
              std::cout << "No pfpart! exit" << std::endl;
              pfpart = false;
              return false;
          }
          
          // and now grab tracks associated to the same PFParts
          larlite::AssSet_t ass_pfpart_trk_v;
          larlite::event_track *ev_trk_2 = nullptr;
          bool tracks = true; // have the tracks been found?
          if (pfpart){
              ass_pfpart_trk_v = storage->find_one_ass( ev_pfpart->id(), ev_trk_2, ev_pfpart->name() );
              
              if (!ev_trk_2 or (ev_trk_2->size() == 0)){
                  std::cout << "No track associated to PFPart! exit" << std::endl;
                  tracks = false;
              }
          }
          
          if (_verbose)
              std::cout << "Associations between vtx and track : " << ass_vtx_trk_v.size() << std::endl;
          
          // find the track and vertex associated to the neutrino
          for (size_t i=0; i < ass_vtx_trk_v.size(); i++){
              
              if (ass_vtx_trk_v[i].size() == 0){
                  std::cout << "vtx->trk association is empty..." << std::endl;
                  continue;
              }
              if (_verbose){
                  std::cout << "trk " << i << " associated to vtx " << ass_vtx_trk_v[i][0] << std::endl;
                  std::cout << ev_trk->size() << " tracks present.." << std::endl;
                  std::cout << ev_vtx->size() << " vertices present.." << std::endl;
              }
              auto const& nutrk = ev_trk->at(i);                    //note this is the muon track (despite the label)
              auto const& nuvtx = ev_vtx->at( ass_vtx_trk_v[i][0] );
              
              std::cout << " nu vtx " << nuvtx.X() << ", " << nuvtx.Y() << ", " << nuvtx.Z() << std::endl;
              TVector3 vtx_pandora;
              vtx_pandora.SetXYZ(nuvtx.X(),nuvtx.Y(),nuvtx.Z());
              double vtx_dist = (vtx_pandora - vtx_pi).Mag();
              std::cout << "vtx_dist" << vtx_dist << std::endl;
              vtx_distance->Fill(vtx_dist);
              double reco_muon_len = nutrk.Length();
              reco_muon_lengths->Fill(reco_muon_len);
              
              // grab the PFParticle associated with this muon
              if (ass_trk_pfpart_v.size() <= i)
                  return false;
              if (_verbose)
                  std::cout << "PFParts associated with muon : " <<  ass_trk_pfpart_v[i].size() << std::endl;
              auto pfpart_idx = ass_trk_pfpart_v[i][0];
              auto muon = ev_pfpart->at(pfpart_idx);
              reco_muontrackstart.SetXYZ(nutrk.LocationAtPoint(0).X(),nutrk.LocationAtPoint(0).Y(), nutrk.LocationAtPoint(0).Z());
              reco_muontrackend.SetXYZ(nutrk.End().X(), nutrk.End().Y(), nutrk.End().Z());
              
              if (_verbose)
                  std::cout << "Muon PFPart info :" << std::endl
                  << "\tPDG code   : " << muon.PdgCode() << std::endl
                  << "\tDaughters? : " << muon.NumDaughters() << std::endl
                  << "\t Parent?   : " << muon.Parent() << std::endl
                  << "\t Length?   : " << nutrk.Length() << std::endl;
              
              // grab parent
              if (muon.Parent() >= ev_pfpart->size()){
                  if (_verbose)
                      std::cout << "Muon parent not here..." << std::endl;
                  return false;
              }
              
              auto neutrino = ev_pfpart->at( muon.Parent() );
              
              if (_verbose)
                  std::cout << "Neutrino PFPart info :" << std::endl
                  << "\tPDG code   : " << neutrino.PdgCode() << std::endl
                  << "\tDaughters? : " << neutrino.NumDaughters() << std::endl
                  << "\t Parent?   : " << neutrino.Parent() << std::endl;
              
              // print neutrino daughters
              for (auto daughter_idx : neutrino.Daughters() ){
                  auto daughter = ev_pfpart->at(daughter_idx);
                  if (_verbose)
                      std::cout << "daughter PFPart info :" << std::endl
                      << "\tPDG code   : " << daughter.PdgCode() << std::endl
                      << "\tDaughters? : " << daughter.NumDaughters() << std::endl
                      << "\t Parent?   : " << daughter.Parent() << std::endl;
                  if (daughter.PdgCode() == 11)
                      n_showers += 1;
                  
                  if (daughter.PdgCode() == 13) {
                      
                      if (tracks == false) // don't search for reco tracks if they have not been found.
                          continue;
                      n_tracks += 1;
                      auto trk_idx = ass_pfpart_trk_v[daughter_idx];
                      if (trk_idx.size() != 1){
                          std::cout << "no associated track...skip" << std::endl;
                          continue;
                      }
                      auto trk = ev_trk_2->at(trk_idx[0]);
                      std::cout << "Added daughter to neutrino tracks" << std::endl;
                      std::cout << "Daughter Track length " << trk.Length() << std::endl;
                      TVector3 trkStart;
                      TVector3 trkEnd;
                      trkStart.SetXYZ(trk.LocationAtPoint(0).X(), trk.LocationAtPoint(0).Y(), trk.LocationAtPoint(0).Z());
                      trkEnd.SetXYZ(trk.End().X(), trk.End().Y(), trk.End().Z());
                      std::cout << "Daughter start Point \t" << trkStart.X() << ", " << trkStart.Y() <<", " << trkStart.Z() << std::endl;
                      std::cout << "Dautghter end Point \t" << trkEnd.X() << ", " << trkEnd.Y() << ", " << trkEnd.Z() << std::endl;
                      double trkvtxdist = -99999 ;
                      if((trkStart - vtx_pandora).Mag() < (trkEnd - vtx_pandora).Mag() ) {
                          trkvtxdist = (trkStart - vtx_pandora).Mag();
                      }
                      else trkvtxdist = (trkEnd - vtx_pandora).Mag();
                      TrackVertexDist->Fill(trkvtxdist);
                      if(trkvtxdist < 25) n_close_track++;
                      // add this track to the list of track-outputs of the neutrino interaction
                      nu_trk_v.push_back( trk );
                      
                  }// if tracks
              }
              
              std::cout << "Number of associated neutrino daughter tracks " << nu_trk_v.size() << std::endl;
              
              
              reco_ntracks->Fill(n_tracks);
              reco_nshowers->Fill(n_showers);
              diff_ntracks->Fill(n_t+n_m - n_tracks);
              diff_nshowers->Fill(n_s_photon - n_showers);
              diff_nclosetracks->Fill(n_t+n_m - n_close_track);
              diff_muon_length->Fill(muon_length - reco_muon_len);
              if(muon_length > 75) diff_muon_length_cut->Fill(muon_length - reco_muon_len);
              if(muon_length < 75) vtx_distance_cut->Fill((vtx_pandora - vtx_pi).Mag());
              dist_muonStart->Fill((muontrackstart - reco_muontrackstart).Mag());
              if(muon_length > 75) dist_muonStart_cut->Fill((muontrackstart - reco_muontrackstart).Mag());
              total_ntrk_nshw->Fill((n_t + n_m + n_s_photon) - (n_close_track + n_showers));
              
              // found the muon -> so exit track loop...
              
              if (_verbose)
                  std::cout << std::endl << std::endl << std::endl;
              
              break;
          }
          
          //if (_filter_showers and (n_showers != _n_showers) )
              return false;
          
          //if (_filter_tracks and (n_tracks != _n_tracks) )
              return false;
          
          return true;
        }
  }
  // ============================================================================================= //

  bool numuCCfilterstudy::finalize() {
      std::cout << "\n =========================================================================== \n " << std::endl;
      std::cout << "\n number of all events processed: \t\t" << all_events << std::endl;
      std::cout << " number of numuCCPi0 events in TPC \t\t" << numuCCpi0_events << std::endl;
      std::cout << " number of noPi0 events             \t\t" << noPi0 << std::endl;
      std::cout << " number of more than one Pi0 events \t\t" << twoPi0<< std::endl;
      std::cout << " muon of events with MC topology    \t\t" << nTopoCut << std::endl;

    // This function is called at the end of event loop.
    // Do all variable finalization you wish to do here.
    // If you need, you can store your ROOT class instance in the output
    // file. You have an access to the output file through "_fout" pointer.
    //
    // Say you made a histogram pointer h1 to store. You can do this:
    //
      if(_fout) { _fout->cd();
          mc_track_lengths->Write();
          mc_muon_lengths->Write();
          mc_ntracks->Write();
          mc_nshowers->Write();
          vtx_distance->Write();
          reco_ntracks->Write();
          reco_nshowers->Write();
          diff_ntracks->Write();
          diff_nshowers->Write();
          reco_muon_lengths->Write();
          diff_muon_length->Write();
          diff_muon_length_cut->Write();
          vtx_distance_cut->Write();
          TrackVertexDist->Write();
          diff_nclosetracks->Write();
          dist_muonStart->Write();
          dist_muonStart_cut->Write();
          total_ntrk_nshw->Write();
          
      }
    //
    // else 
    //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
    //
  
    return true;
  }

}
#endif
