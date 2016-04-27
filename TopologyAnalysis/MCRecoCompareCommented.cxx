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
      
      /*
      if(p_tracks_v->size() == 0) std::cout << "\n no reco tracks found" << std::endl;
      else std::cout << "\n number of reco tracks " << p_tracks_v->size() << std::endl;
      
      if(p_showers_v->size() == 0) std::cout << " no reco showers found" << std::endl;
      else std::cout << "\n number of reco showers " << p_showers_v->size() << std::endl;
      */
      
      
      // Track associations - - - - - - - - - - - - - - - - - -
      
      auto *pfp_v = storage->get_data<event_pfpart>("pandoraNu");               // get PFParticles
      /*
      event_track *pfp_tracks = nullptr;
      auto ass_tracks = storage->find_one_ass(pfp_v->id(), pfp_tracks, pfp_v->name() ); // this *should* be the tracks associated with the pfparticle...
      
      std::cout << "ass_tracks size " << ass_tracks.size() <<std::endl;
      */
      //
      
      
      int pTrk      = 0;        // number of PFParticles that are Tracks
      int pShow     = 0;        // number of PFParticles that are Showers
      bool p_topoCut = false;    // flag for 2 track, 2 shower events (initial seeding out)
      
      
      for(auto pfp : *pfp_v) {
          if(pfp.PdgCode() == 12 || pfp.PdgCode() == 14) {
              //std::cout <<"\n Neutrino index: " << pfp.Self() << std::endl;
          }
          
          if(abs(pfp.PdgCode()) == 13) {
              pTrk++;
              //std::cout << " track pfp.Parent() index " << pfp.Parent() << "\t parent pdg " << pfp_v->at(pfp.Parent()).PdgCode()<< std::endl;
          }
          else if(pfp.PdgCode() == 11) {
              pShow++;
              //std::cout << " shower pfp.Parent() index " << pfp.Parent() << "\t parent pdg " << pfp_v->at(pfp.Parent()).PdgCode()<< std::endl;
          }
          //else std::cout << "neither shower nor track " << pfp.PdgCode() << std::endl;
      }
      
      
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
          }
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
      
      int nMuTrk    = 0;        //number of primary muon tracks
      int nTrk      = 0;        //number of tracks other than primary muons (may include non-primary muons)
      int prim_trk  = 0;        //number of primary tracks other than muons (what else comes out of int)
      int nShow     = 0;        //number of showers
      int prim_show = 0;        //number of showers from the pi0 decay
      
      
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
      if(nPi0 >= 1 && nMu == 1) {
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
              if(mcp.PdgCode() == 111 && mcp.StatusCode() == 1) {                                                                                // find pizero
                  std::cout << "FOUND PIZERO " << mcp.TrackId()  << " pos " << mcp.Position().X() << ",  \t daughters: " << mcp.Daughters().size() << std::endl;
                  int pizID = mcp.TrackId();
                  for(auto mcshw : *mcshow_v) {
                      //if(mcshw.PdgCode() == 22) std::cout << "\tPhoton shower " << mcshw.MotherTrackID() << ", " << mcshw.AncestorTrackID() << " starts " << mcshw.Start().X() << std::endl;
                      if(mcshw.PdgCode() == 22) {
                          n_s++;           //count photon showers from pizero
                          std::cout << "\tPhoton shower " << mcshw.MotherTrackID() << ", " << mcshw.AncestorTrackID() << " starts " << mcshw.Start().X() << std::endl;
                      }
                  }
                  for(auto mctrk : *mctrk_v) {
                      if(abs(mctrk.PdgCode()) == 13 && mctrk.Process() == "primary") {
                          n_m++;         // primary muons
                          std::cout << "\t Muon " << std::endl;
                      }
                      else if(mctrk.PdgCode() != 2112 && mctrk.PdgCode() != 130 && mctrk.PdgCode() != 310 && mctrk.PdgCode() != 111
                              && mctrk.PdgCode() != 311 && mctrk.PdgCode() < 10000 && mctrk.Process() == "primary") {
                                n_t++;           // other tracks from interaction
                                std::cout << "\t Other " << mctrk.PdgCode() << std::endl;
                      }
                  }
                  std::cout << "\t n_s " << n_s << "\t n_m " << n_m << "\t n_t " << n_t << std::endl;
                  if(n_s == 2 && n_m == 1 && n_t >=1) {                                                                             // selection Criteria
                      topoCut = true;
                      nTopoCut++;
                  }
              }
          }
          
          
      } // selecting numuCCpi0 events
      
      
     /*
          
          for(auto mctrk : *mctrk_v){
              if(StartInTPC(mctrk) == true){                                                                        //ensure looking at tracks inside TPC volume
                  
                  if(abs(mctrk.PdgCode()) == 13 && mctrk.MotherPdgCode() == 111 && mctrk.MotherTrackID() == mctrk.AncestorTrackID()) {    // primary muon - out of PiZero interaction
                      TVector3 start  = { mctrk.Start().X(), mctrk.Start().Y(), mctrk.Start().Z()};
                      TVector3 end    = { mctrk.End().X(),mctrk.End().Y(),mctrk.End().Z() };
                      muon_length = length(mctrk);                                                                      // using length function - mcsteps added up
                      vtx = start;                                                                                      // vtx: start of primary muon,from numu interaction
                      //std::cout << " \n MC vertex at " << vtx.X() << ","<<vtx.Y()<<","<<vtx.Z()<<std::endl;
                      //std::cout << " MC muon length " << muon_length << std::endl;
                      int n_s = 0;
                      int n_t = 0;
                      for(mctrkk : *mctrk_v) {
                          if(mctrkk.PdgCode() != 2112 && mctrkk.PdgCode() != 111 & mctrkk.PdgCode() != 130 && mctrkk.PdgCode() != 310 &&
                             mctrkk.PdgCode() != 311 & mctrkk.PdgCode() < 10000 && mctrkk.MotherTrackID() == mctrk.MotherTrackID) {
                              n_t++;
                          }
                      }
                      for(mcshw : *mcshow_v) {
                          if(mcshw.PdgCode() == 22 && mcshw.MotherPdgCode() == 111) n_s++;
                      }
                      if(n_s == 2 && n_t >= 1) topoCut = true;
                  }
                  
           
                  }
                  else if(mctrk.PdgCode() != 2112 && mctrk.PdgCode() != 111 && mctrk.PdgCode() != 130
                          && mctrk.PdgCode() != 310 && mctrk.PdgCode() != 311 && mctrk.PdgCode() < 10000){            // cutting out neutral particles and fragments
                      TVector3 start  = { mctrk.Start().X(), mctrk.Start().Y(), mctrk.Start().Z()};
                      TVector3 end    = { mctrk.End().X(),mctrk.End().Y(),mctrk.End().Z() };
                      other_length = length(mctrk);
                      d_vtx_trk = (vtx - start).Mag(); // distance to vertex
                      other_pdg = mctrk.PdgCode(); //std::cout << other_pdg << std::endl;
                      other_process = mctrk.Process();
                      if(other_process == "primary" && other_length != 0) prim_trk++;                                   // making sure track is not zero length
                     // std::cout << "\n other MC track in event " << mctrk.Start().X()<<", "<< mctrk.Start().Y()<<", " << mctrk.Start().Z()<< std::endl;
                     // std::cout << " other MC track length " << other_length << std::endl;
                      if(other_length != 0) nTrk++;                                                                      // this counts tracks other than muons
                  }
                
            }
          } //acting on each element of mctrk_v
          
           
           
          for(auto mcshow : *mcshow_v) {
              if(StartInTPC(mcshow) == true) {
                  nShow++;
                  if(mcshow.MotherPdgCode() == 111 && mcshow.MotherTrackID() == mcshow.AncestorTrackID()) {
                      prim_show++;
                      TVector3 start  = { mcshow.Start().X(), mcshow.Start().Y(), mcshow.Start().Z()};
                      d_vtx_show = (vtx - start).Mag();
                     // std::cout << " \n MC shower from " <<mcshow.Start().X()<<","<<mcshow.Start().Y()<<","<<mcshow.Start().Z()<< std::endl;
                  }
              }
          } // acting on each element of mcshow_v
          
          // count events that pass the topoCut
          if(prim_trk >= 1 && nMuTrk == 1 && prim_show == 2) {
              topoCut = true;                                                   // topoCut is 1 muon track, AT LEAST one other track and 2 showers from vertex
              nTopoCut++;
          }
          */
          
          
      
           
    
      
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
      
      if(p_topoCut && topoCut == false) {
          nFake++;
      }
      
      /*    std::cout << "\n Event: " << storage->event_id() << std::endl;
          for(auto pfp : *pfp_v) {
              std::cout << " PFParticles in fake events " << pfp.PdgCode() << " parent : " << pfp.Parent()<<std::endl;
          }
          std::cout << "\n" << std::endl;
          for(auto mctrk : *mctrk_v){
              std::cout << " MCTracks in fake events " << mctrk.PdgCode() << " \t own track ID " << mctrk.TrackID() << "\t Mother PDG " << mctrk.MotherPdgCode() << "\t Ancestor and Mother ID " << mctrk.MotherTrackID() << ", " <<mctrk.AncestorTrackID();
              if(StartInTPC(mctrk) == true)std::cout << " \t StartInTPC true ";
              std::cout << std::endl;
          }
          for(auto mcshw : *mcshow_v) {
              std::cout << " MCShowers in fake events " << mcshw.PdgCode() << "\t own track ID " << mcshw.TrackID() << "\t Mother PDG " << mcshw.MotherPdgCode() << "\t Ancestor and Mother ID " << mcshw.MotherTrackID() << ", " <<mcshw.AncestorTrackID();
              if(StartInTPC(mcshw) == true)std::cout << " \t StartInTPC true ";
              std::cout << std::endl;
          }
      } */
  
    return true;
  }
    
//=================================================================================================//
//=================================================================================================//

  bool MCRecoCompare::finalize() {
      
      std::cout << "\n number of all events processed: \t\t" << all_events << std::endl;
      std::cout << " number of CC events \t\t" << CC_events << std::endl;
      std::cout << " number of numuCCPi0 events \t\t" << all << std::endl;
      std::cout << " number passing MC topoCut  \t\t" << nTopoCut << std::endl;
      std::cout << " number passing p_topoCut \t\t" << p_nTopoCut << std::endl;
      
      std::cout << "topoCut and p_topoCut true \t\t" << nSuccess << std::endl;
      std::cout << "topoCut false and p_topoCut true \t\t" << nFake << std::endl;
    
      
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
