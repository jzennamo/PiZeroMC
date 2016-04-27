#ifndef LARLITE_MCRECOCOMPARE_CXX
#define LARLITE_MCRECOCOMPARE_CXX

#include "MCRecoCompare.h"

namespace larlite {

  bool MCRecoCompare::initialize() {
  
  // General
      all_events      = 0;
      
      
  // FOR MC
      
      
      evt_id = 0; ///only for debug
      
      CC_events         = 0;
      all               = 0; //these are all CCpi0 events with one muon, one pi0 (numuCCPi0)
      nTopoCut          = 0; //number of events passing topoCut
      
      if(tree) delete tree;
      tree = new TTree("tree","tree");
      tree->Branch("thisInt",&thisInt,"thisInt/I");
      tree->Branch("lengthRatio", &lengthRatio, "lengthRatio/D");   // Ratio of other track length to muon track length
      
      muons = new TTree("muons","muons");
      muons->Branch("muon_length",&muon_length,"muon_length/D");    // length of muon in event that passes the topoCut
      
      other = new TTree("other","other");
      other->Branch("other_length",&other_length,"other_length/D"); // length of other tracks in event that passes topoCut
      
      fake = new TTree("fake","fake");
      
      

  // FOR RECO
      
      p_CCpi0_events     = 0;
      p_all              = 0;
      p_nTopoCut         = 0;
      
      nSuccess           = 0; // number of MC and reco matches
      nFake              = 0;
      nCCpi0Pand     = 0; // number of events with p_topoCut true and numuCCpi0 true

    return true;
  }
  
//=================================================================================================//
//=================================================================================================//
    
  bool MCRecoCompare::analyze(storage_manager* storage) {
      

      //std::cout << "\n ======================================================== " <<std::endl;
      //std::cout << " EVENT NUMBER " << storage->event_id() << std::endl;
      
      all_events++;
      
      // Reco data products ------------------------------------------------------------------------------
      
      //auto p_tracks_v = storage->get_data<event_track>("pandoraNuKHit");        // get reco tracks
      //auto p_showers_v = storage->get_data<event_shower>("showerrecopandora");  // get reco showers
      
      
      // Track associations - - - - - - - - - - - - - - - - - -
      
      auto *pfp_v = storage->get_data<event_pfpart>("pandoraNu");               // get PFParticles
      
      event_track *pfp_tracks = nullptr;
      auto ass_tracks = storage->find_one_ass(pfp_v->id(), pfp_tracks, pfp_v->name() ); // this *should* be the tracks associated with the pfparticle...
      
      //std::cout << "ass_tracks size " << ass_tracks.size() <<std::endl;
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
      
      bool p_topoCut = false;    // flag for 2 track, 2 shower events (initial seeding out)
      
      
      
      // go through PFParticle list, check for neutrino (12 or 14)
      for(auto pfp : *pfp_v) {
          if(pfp.PdgCode() == 12 || pfp.PdgCode() == 14) {                                      //is PFParticle a neutrino?
              //std::cout << "\n neutrino " << pfp.PdgCode();
              int nu_index = pfp.Self(); //std::cout << nu_index << std::endl;                    // if yes, this is its identifier
              int n_t = 0;                                                                          // number of tracks
              int n_s = 0;                                                                          // number of showers
              for(auto pfpp : *pfp_v) {
                  if(pfpp.PdgCode() == 13 && pfpp.Parent()== nu_index) n_t++;                          // number of tracks from that neutrino
                  if(pfpp.PdgCode() == 11 && pfpp.Parent()== nu_index) n_s++;                          // number of showers from that neutrino
              }
            // std::cout << "\n" << n_s <<" SHOWERS AND " << n_t<< " TRACKS \n\n\n" << std::endl;
              if(n_s == 2 && n_t >= 2) {
                  //std::cout << "\n 2 TWO SHOWERS AND 2 TRACKS!! \n" << std::endl;
                  p_nTopoCut++;
                  p_topoCut = true;
              }
              //else std::cout << "\t Not right number of showers and tracks" <<std::endl;
          }
          //else std::cout << "No neutrino found " <<std::endl;
      }
      
      
      //MCTruth data products -----------------------------------------------------------------------------
      auto mc = storage->get_data<event_mctruth>("generator");
      auto& all_particles = mc->at(0).GetParticles();
      auto& nu = mc->at(0).GetNeutrino();
      
      int nMu       = 0;        //number of MCTruth muons
      int nPi0      = 0;        // number of MCTruth Pi0
      
      bool numuCCpi0 = false;   //flag for numuCCpi0 event
      bool topoCut   = false;   //flag for topology cut
      
      
      //MCTrack / MCShower data products ------------------------------------------------------------------
      
      
      //MCTracks and MCShowers
      auto mctrk_v = storage->get_data<event_mctrack>("mcreco");        // this is a vector of tracks
      auto mcshow_v = storage->get_data<event_mcshower>("mcreco");      // this is a vector of showers
      
      
      // Check whether event is CC numu
      if(nu.CCNC() == 0 && abs(nu.Nu().PdgCode()) == 14) {
          CC_events++;
          for(auto const& par: all_particles) {
              if(par.StatusCode() == 1) {                                                                           //look at stable final states only
                  if(par.PdgCode() == 111)      nPi0++;
                  if(abs(par.PdgCode()) == 13)  nMu++;
              }
          }
      } // select numu CC event
      
      
      //these are the events we want to initially select and check their tracks and showers, numuCCPi0
      if(nPi0 >= 1 && nMu >= 1) {
          numuCCpi0 = true;
          all++;
      }
      
      //////////////////////////////////////////////////////////////////////////////////
      // This bit only executed when looking at a numuCCpi0 event                     //
      
     // std::cout << "\n event " << std::endl;
      
      if(numuCCpi0) {
          
          
          // Vertex
          TVector3 vtx = {-9999,-9999,-9999};
          
          for(auto mcp : all_particles) {                                                                                      //iterate through MCParticles
              int n_s = 0;
              int n_m = 0;
              int n_t = 0;
              double muon_len = 0;
              if(mcp.PdgCode() == 111 && mcp.StatusCode() == 1) {                                                                                // find pizero
                  //std::cout << "FOUND PIZER0, mother:  " << mcp.Mother() << std::endl;
                  //int pizero_mother = mcp.Mother();
                  for(auto mcshw : *mcshow_v) {
                      //if(mcshw.PdgCode() == 22) std::cout << "\tPhoton shower " << mcshw.MotherTrackID() << ", " << mcshw.AncestorTrackID() << " starts " << mcshw.Start().X() << std::endl;
                      if(mcshw.PdgCode() == 22) {
                          n_s++;           //count photon showers from pizero
                          //std::cout << "\tPhoton shower " << mcshw.MotherTrackID() << ", " << mcshw.AncestorTrackID() << " starts " << mcshw.Start().X() << std::endl;
                      }
                  }
                  for(auto mctrk : *mctrk_v) {
                      if(abs(mctrk.PdgCode()) == 13 && mctrk.Process() == "primary") {
                          if(length(mctrk) > 0) { n_m++;         // primary muons
                              muon_len = length(mctrk);}
                          //std::cout << "\t Muon " << std::endl;
                      }
                      else if(mctrk.PdgCode() != 2112 && mctrk.PdgCode() != 130 && mctrk.PdgCode() != 310 && mctrk.PdgCode() != 111
                              && mctrk.PdgCode() != 311 && mctrk.PdgCode() < 10000 && mctrk.Process() == "primary") {
                                n_t++;           // other tracks from interaction
                                //std::cout << "\t Other " << mctrk.PdgCode() << std::endl;
                      }
                  }
                  //std::cout << "\t n_s " << n_s << "\t n_m " << n_m << "\t n_t " << n_t << std::endl;
                  if(n_s == 2 && n_m >= 1 && n_t >=1) {                                                                             // selection Criteria
                      topoCut = true;
                      nTopoCut++;
                      for(auto mck : all_particles) {
                          //f(mck.TrackId() == pizero_mother) std::cout << " pizero mother " << mck.PdgCode() << "\t muon length " << muon_len<< std::endl;
                      }
                  }
                  else if(p_topoCut) std::cout << "p_topoCUt true but not right number of MCTracks \t muons: "<< n_m << "\t tracks: " << n_t << "\t showers: " << n_s << std::endl;
              }
            
          }
          
          
      } // selecting numuCCpi0 events
      

    
      
      //////////////////////////////////////////////////////////////////////////////////////////
    /*
      lengthRatio = other_length/muon_length;
      
      //Filling appropriate trees ==== NEEDS WORK
      if(topoCut) {
          muons->Fill();
          other->Fill();
          tree->Fill();
      }

      */
    
     if(topoCut && p_topoCut) {
          nSuccess++;
          //std::cout << storage->event_id() << std::endl;
      }
      
      if(p_topoCut && numuCCpi0) nCCpi0Pand++;
      
      if(p_topoCut && topoCut == false) nFake++;
      
      if(p_topoCut && topoCut == false) {
          //std::cout << "neutrino " << nu.Nu().PdgCode() <<std::endl;
          for(auto mcp : all_particles) {}//std::cout << " mc particles " << mcp.PdgCode() << std::endl;
      }
      
    return true;
  }
    
//=================================================================================================//
//=================================================================================================//

  bool MCRecoCompare::finalize() {
      
      std::cout << "\n number of all events processed: \t\t" << all_events << std::endl;
      std::cout << " number of CC events \t\t\t\t" << CC_events << std::endl;
      std::cout << " number of numuCCPi0 events \t\t\t" << all << std::endl;
      std::cout << " number passing MC topoCut  \t\t\t" << nTopoCut << std::endl;
      std::cout << " number passing p_topoCut \t\t\t" << p_nTopoCut << std::endl;
      
      std::cout << " MC topoCut true and p_topoCut true \t\t" << nSuccess << std::endl;
      std::cout << " MC topoCut false and p_topoCut true \t\t" << nFake << std::endl;
      std::cout << "\n  numuCCPi0 event and p_topoCut true \t\t" << nCCpi0Pand << std::endl;
    
      
      // This function is called at the end of event loop.
      // Do all variable finalization you wish to do here.
      // If you need, you can store your ROOT class instance in the output
      // file. You have an access to the output file through "_fout" pointer.
      //
      // Say you made a histogram pointer h1 to store. You can do this:
      //
      if(_fout) { _fout->cd(); tree->Write(); muons->Write(); other->Write();}
      
      //else
      //   print(MSG::ERROR,__FUNCTION__,"Did not find an output file pointer!!! File not opened?");
      //
      
    return true;
  }

}
#endif
