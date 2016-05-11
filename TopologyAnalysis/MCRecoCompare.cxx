#ifndef LARLITE_MCRECOCOMPARE_CXX
#define LARLITE_MCRECOCOMPARE_CXX

#include "MCRecoCompare.h"

namespace larlite {

  bool MCRecoCompare::initialize() {
  
  // General
      all_events      = 0;                  // count all processed events
      
      
  // FOR MC
      
      
      evt_id = 0; ///only for debug
      
      CC_events         = 0;
      numuCCpi0_events  = 0; //these are all CCpi0 events with one muon, one pi0 (numuCCPi0)
      
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
      
      p_nTopoCut         = 0;       // number of reco events passing p_topoCut
      
  // Evaluation
      
      nSuccess           = 0;       // number of MC and reco matches
      nFake              = 0;       // number of reco PiZero events where there was no MC PiZero event
      fakeNotCC          = 0;       // number of fakes (true event not CC)
      fakeNotTopo        = 0;       // number of fakes (MC doesnt have the required topology (number of showers, tracks))
      

    return true;
  }
  
//=================================================================================================//
//=================================================================================================//
    
  bool MCRecoCompare::analyze(storage_manager* storage) {
      

      //std::cout << " ======================================================== " <<std::endl;
      //std::cout << " EVENT NUMBER " << storage->event_id() << std::endl;
      
      all_events++;

      
      // Reco data products ------------------------------------------------------------------------------
      
      //auto *p_tracks_v = storage->get_data<event_track>("pandoraNu");                // get reco tracks
      //auto p_showers_v = storage->get_data<event_shower>("showerrecopandora");    // get reco showers
      
      
      // Track associations - - - - - - - - - - - - - - - - - -
      
      auto *pfp_v = storage->get_data<event_pfpart>("pandoraNu");               // get PFParticles
      
      //event_track *pfp_tracks = nullptr;
      //auto ass_tracks = storage->find_one_ass(pfp_v->id(), pfp_tracks, pfp_v->name() ); // this *should* be the tracks associated with the pfparticle...
      
      //std::cout << "ass_tracks size " << ass_tracks.size() <<std::endl;
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
      bool moreThanOneNu = false;   // flag event if there are more than one reco neutrinos in the event
      bool foundOne      = false;   // internal flag
      bool p_topoCut     = false;   // flag for 2 track, 2 shower events (initial seeding out)
      
      
      // go through PFParticle list, check for neutrino (12 or 14)
      
      for(auto pfp : *pfp_v) {
          if((pfp.PdgCode() == 12 || pfp.PdgCode() == 14) && foundOne == false) {                    //is PFParticle a neutrino? Only go on if not already have a success
              int nu_index = pfp.Self(); //std::cout << nu_index << std::endl;                       // if yes, this is its identifier
              int n_t = 0;                                                                           // number of tracks
              int n_s = 0;                                                                           // number of showers
              for(auto pfpp : *pfp_v) {
                  if(pfpp.PdgCode() == 13 && pfpp.Parent()== nu_index) n_t++;                          // number of tracks from that neutrino
                  if(pfpp.PdgCode() == 11 && pfpp.Parent()== nu_index) n_s++;                          // number of showers from that neutrino
              }
              if(n_s == 2 && n_t >= 2) {
                  p_topoCut = true;                                                                     // p_topoCut is flag for reco event -> query to count at the end.
                  foundOne = true;
              }
          }
          else if((pfp.PdgCode() == 12 || pfp.PdgCode() == 14) && foundOne == true) {               // if already have a selected reco event
              int nu_index = pfp.Self(); //std::cout << nu_index << std::endl;
              int n_t = 0;
              int n_s = 0;
              for(auto pfpp : *pfp_v) {
                  if(pfpp.PdgCode() == 13 && pfpp.Parent()== nu_index) n_t++;
                  if(pfpp.PdgCode() == 11 && pfpp.Parent()== nu_index) n_s++;
              }
              if(n_s == 2 && n_t >= 2) moreThanOneNu = true;                                        // more than one reco event.
          }
      }

      if(p_topoCut && moreThanOneNu == false) p_nTopoCut++;
      
      
      
      
      
      //MCTruth data products -----------------------------------------------------------------------------
      auto mc = storage->get_data<event_mctruth>("generator");      // these are MCParticles from generator level
      auto& all_particles = mc->at(0).GetParticles();
      auto& nu = mc->at(0).GetNeutrino();
      
      int nMu       = 0;                //number of MCTruth muons
      int nPi0      = 0;                // number of MCTruth Pi0
      
      bool numuCCpi0        = false;    //flag for numuCCpi0 event
      bool topoCut          = false;    //flag for topology cut
      
      
      //MCTrack / MCShower data products ------------------------------------------------------------------
      
      //MCTracks and MCShowers
      auto mctrk_v = storage->get_data<event_mctrack>("mcreco");        // this is a vector of MCtracks
      auto mcshow_v = storage->get_data<event_mcshower>("mcreco");      // this is a vector of MCshowers
      
      //auto ev_mcpart = storage->get_data<event_mcpart>("mcreco");     // this would be MCParticle from geant, but not in event
      //event_mctrack * ev_mctrack = nullptr;
      //auto track_assoc = storage->find_one_ass(ev_mcpart->id(), ev_mctrack, ev_mcpart->name());
      
      //if(mc->at(0).Origin() != 1) std::cout << "origin " << mc->at(0).Origin() << std::endl;
    
      // Check whether event is CC numu
      
      TVector3 vtx = {-9999,-9999,-9999}; // Pion Vertex

      
      if(nu.CCNC() == 0 && (abs(nu.Nu().PdgCode()) == 14) ) {
          CC_events++;
          for(auto const& par: all_particles) {
              if(par.PdgCode() == 111 && StartInTPC(par) && par.StatusCode() == 1) {                //Pions in TPC only
                nPi0++;
                vtx.SetXYZ(par.Trajectory().at(0).Position().X(),par.Trajectory().at(0).Position().Y(),par.Trajectory().at(0).Position().Z());
              }
          }
          if(nPi0 == 1) {
              for(auto const& par : all_particles) {
                  if(abs(par.PdgCode()) == 13 && StartInTPC(par) && par.Mother() == 0)  {
                      TVector3 muonvtx;
                      muonvtx.SetXYZ(par.Trajectory().at(0).Position().X(),par.Trajectory().at(0).Position().Y(),par.Trajectory().at(0).Position().Z());
                      if((muonvtx - vtx).Mag() < 0.01) nMu++;
                  }
              }
              
          }
          
          if(nPi0 == 1 && nMu == 1) {
              numuCCpi0 = true;
              numuCCpi0_events++;
          }
      } // select numu CC event with one Pion only

      
      
      // now only looking at numuCCpi0 flagged events
      
      if(numuCCpi0) {               // now already know this event has at one Pi0, and one muon coming out the vertex
          
          // Vertex is vtx
          
          std::cout << "\n\n\n\n****** Pion Vertex " << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z() << std::endl;
          
          
          for(auto mcp : all_particles) {                                                                                                       //iterate through MCParticles
          
              if(mcp.PdgCode() == 111 && mcp.StatusCode() == 1 ){                                                                               // find pizero
                  
                  std::cout << storage->event_id()<< "\t" << mcp.PdgCode() << "\t "<< mcp.Trajectory().at(0).Position().X() << ", " << mcp.Trajectory().at(0).Position().Y() << ", " << mcp.Trajectory().at(0).Position().Z() << "\t" << mcp.StatusCode() << "\t Mother ID " << mcp.Mother() << std::endl;
              }
          }
          
          for(auto mcp : all_particles) {
              
              if(abs(mcp.PdgCode()) == 13 && mcp.StatusCode() == 1){
                  std::cout << storage->event_id()<< "\t" << mcp.PdgCode() << "\t "<< mcp.Trajectory().at(0).Position().X() << ", " << mcp.Trajectory().at(0).Position().Y() << ", " << mcp.Trajectory().at(0).Position().Z() << "\t \t Mother ID " << mcp.Mother() << std::endl;
                  TVector3 muonvtx;
                  muonvtx.SetXYZ(mcp.Trajectory().at(0).Position().X(),mcp.Trajectory().at(0).Position().Y(),mcp.Trajectory().at(0).Position().Z());
                  //std::cout << "\t distance to pion vertex " << (vtx - muonvtx).Mag() << std::endl;
              }
              
              
              if(abs(mcp.PdgCode()) != 13 && mcp.PdgCode() != 111 && mcp.PdgCode() != 2112 && mcp.StatusCode() == 1){
                  TVector3 othervtx;
                  othervtx.SetXYZ(mcp.Trajectory().at(0).Position().X(),mcp.Trajectory().at(0).Position().Y(),mcp.Trajectory().at(0).Position().Z());
                  if((vtx - othervtx).Mag() < 10) {
                      std::cout << storage->event_id()<< "\t" << mcp.PdgCode() << "\t "<< mcp.Trajectory().at(0).Position().X() << ", " << mcp.Trajectory().at(0).Position().Y() << ", " << mcp.Trajectory().at(0).Position().Z() << std::endl;
                      //std::cout << "\t distance to pion vertex " << (vtx - othervtx).Mag() << std::endl;
                  }
              }
          }
          
              
              
        // now assess what event looks like in detector: use MCTracks and MCShowers
          std::cout << "\n MCTrack and MCShower profile" << std::endl;
          
            int n_s         = 0;        // number of showers not from Pi0 in vicinity of vertex
            int n_s_photon  = 0;        // number of photon showers from Pi0
            int n_m         = 0;        // number of muons from vertex
            int n_t         = 0;        //number of tracks from vertex
            int n_t_other   = 0;        //number of tracks not from vertex
            double muon_len = 0;
                          
                
          //count tracks
            for(auto mctrk : *mctrk_v) {
                
                if(abs(mctrk.PdgCode()) == 13 && mctrk.Process() == "primary") {
                    TVector3 muontrackstart;
                    muontrackstart.SetXYZ(mctrk.Start().X(), mctrk.Start().Y(), mctrk.Start().Z());
                    if(length(mctrk) > 0 && (muontrackstart - vtx).Mag() < 0.1) {                       // can impose length cut here
                        n_m++;                                                                                            // primary muons
                        muon_len = length(mctrk);
                        std::cout << "*MCTrack muon distance to vtx " << (muontrackstart - vtx).Mag() << "\t len: " << muon_len << std::endl;
                    }
                }
                if(abs(mctrk.PdgCode()) != 13 && mctrk.PdgCode() != 2112 && mctrk.PdgCode() != 130 && mctrk.PdgCode() != 310 && mctrk.PdgCode() != 111
                                      && mctrk.PdgCode() != 311 && mctrk.PdgCode() < 10000 && mctrk.Process() == "primary"  && StartInTPC(mctrk)) {
                    TVector3 trackstart;
                    trackstart.SetXYZ(mctrk.Start().X(), mctrk.Start().Y(), mctrk.Start().Z());
                    if((trackstart - vtx).Mag() < 0.1 && length(mctrk) > 0) {                                                // count other tracks from vertex
                        n_t++;
                        std::cout << "*MCTrack " << mctrk.PdgCode() << " distance to vtx " << (trackstart - vtx).Mag() << "\t len: " << length(mctrk) << std::endl;
                    }
                    // other tracks from interaction
                }
                else if(mctrk.PdgCode() != 2112 && mctrk.PdgCode() != 130 && mctrk.PdgCode() != 310 && mctrk.PdgCode() != 111
                        && mctrk.PdgCode() != 311 && mctrk.PdgCode() < 10000 && mctrk.Process() == "primary"  && StartInTPC(mctrk)) {
                    TVector3 trackstart;
                    trackstart.SetXYZ(mctrk.Start().X(), mctrk.Start().Y(), mctrk.Start().Z());
                    if((trackstart - vtx).Mag() > 0.01 && (trackstart - vtx).Mag() < 5 && length(mctrk) > 0) {                                                // count other tracks NOT from vertex but close
                        n_t_other++;
                        std::cout << "MCTrack " << mctrk.PdgCode() << " distance to vtx " << (trackstart - vtx).Mag() << "\t len: " << length(mctrk) << std::endl;
                    }

                }
            }
          
          
                          
            //count showers
            for(auto mcshw : *mcshow_v) {
                if(mcshw.PdgCode() == 22 && StartInTPC(mcshw) && mcshw.MotherPdgCode() == 111) {                  //count photon showers from pizero
                    double shower_len = 0;
                    TVector3 svtx;
                    TVector3 send;
                    svtx.SetXYZ(mcshw.Start().X(), mcshw.Start().Y(), mcshw.Start().Z());
                    send.SetXYZ(mcshw.End().X(), mcshw.End().Y(), mcshw.End().Z());
                    shower_len = (svtx - send).Mag();
                    std::cout << "**MCShower Photon from Pi0 distance to vtx : " << (svtx - vtx).Mag() << "\t shower length " << shower_len << "\t Mother Process: " << mcshw.MotherProcess() << std::endl;
                    if(shower_len > 0) n_s_photon++;                   //can impose shower length cut here
                }
                else if((mcshw.PdgCode() == 22 || abs(mcshw.PdgCode()) == 11) && StartInTPC(mcshw)) {             // count other showers
                    TVector3 svtx;
                    svtx.SetXYZ(mcshw.Start().X(), mcshw.Start().Y(), mcshw.Start().Z());
                    if((svtx - vtx).Mag() < 100){
                        n_s++;                                                             // count other showers that are close to vertex
                        std::cout << "MCShower " << mcshw.PdgCode() << " distance to vtx : " << (svtx - vtx).Mag() << "\t Mother Process " << mcshw.MotherProcess() << std::endl;
                    }
                }
            }
                          
                          
                          // assess topology
                          //std::cout << "\t n_s " << n_s << "\t n_m " << n_m << "\t n_t " << n_t << std::endl;
                          
                if(n_s_photon == 2 && n_m == 1 && n_t >=1 && n_s == 0 && n_t_other == 0) {                                                                    // selection Criteria
                    topoCut = true;
                            
                          }
                  
              
          
        
    } // selecting numuCCpi0 events

      if(p_topoCut && numuCCpi0 == false && topoCut == false) fakeNotCC++;
      if(p_topoCut && numuCCpi0 && topoCut == false) fakeNotTopo++;
    
      if(topoCut == true) nTopoCut++;

      
      if(topoCut == true && p_topoCut == true) nSuccess++;
      if(p_topoCut && topoCut == false) nFake++;
      
      
      if(numuCCpi0) std::cout << "\n === topoCut : " << topoCut << "\t \t p_topoCut : " << p_topoCut << std::endl;
      
    return true;
  }
    
//=================================================================================================//
//=================================================================================================//

  bool MCRecoCompare::finalize() {
      
      std::cout << "\n number of all events processed: \t\t" << all_events << std::endl;
      std::cout << " number of CC events \t\t\t\t" << CC_events << std::endl;
      std::cout << " number of numuCCPi0 events in TPC \t\t" << numuCCpi0_events << std::endl;
      std::cout << " number passing MC topoCut  \t\t\t" << nTopoCut << std::endl;
      std::cout << " number passing p_topoCut \t\t\t" << p_nTopoCut << std::endl;
      
      std::cout << " MC topoCut true and p_topoCut true \t\t" << nSuccess << std::endl;
      std::cout << " MC topoCut false and p_topoCut true (Fake) \t" << nFake << std::endl;
      std::cout << " \t of which not numuCCpi0 \t\t" << fakeNotCC << std::endl;
      std::cout << " \t of which not correct topology in MC \t" << fakeNotTopo << std::endl;
    
    
      
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
