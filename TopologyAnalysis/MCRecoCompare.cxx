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
      nMoreThanOne       = 0;
      

    return true;
  }
  
//=================================================================================================//
//=================================================================================================//
    
  bool MCRecoCompare::analyze(storage_manager* storage) {
      

      //std::cout << " ======================================================== " <<std::endl;
      //std::cout << " EVENT NUMBER " << storage->event_id() << std::endl;
      
      all_events++;
      
      // Reco data products ------------------------------------------------------------------------------
      
      auto *p_tracks_v = storage->get_data<event_track>("pandoraNu");                // get reco tracks
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
      auto mc = storage->get_data<event_mctruth>("generator");
      auto& all_particles = mc->at(0).GetParticles();
      auto& nu = mc->at(0).GetNeutrino();
      
      int nMu       = 0;                //number of MCTruth muons
      int nPi0      = 0;                // number of MCTruth Pi0
      
      bool numuCCpi0        = false;    //flag for numuCCpi0 event
      bool topoCut          = false;    //flag for topology cut
      bool foundOneMC       = false;    //internal flag
      bool moreThanOneNuMC  = false;    //flag event if there are more than one MC neutrino PiZero interactions
      
      
      //MCTrack / MCShower data products ------------------------------------------------------------------
      
      //MCTracks and MCShowers
      auto mctrk_v = storage->get_data<event_mctrack>("mcreco");        // this is a vector of tracks
      auto mcshow_v = storage->get_data<event_mcshower>("mcreco");      // this is a vector of showers
      
      
      // Check whether event is CC numu
      if(nu.CCNC() == 0 && (abs(nu.Nu().PdgCode()) == 14) ) {
          CC_events++;
          for(auto const& par: all_particles) {
              if(par.StatusCode() == 1) {                                //look at stable final states only
                  if(par.PdgCode() == 111)      nPi0++;
                  if(abs(par.PdgCode()) == 13)  nMu++;
              }
          }
          
          if(nPi0 >= 1 && nMu >= 1) {
              numuCCpi0 = true;
              numuCCpi0_events++;
          }
      } // select numu CC event with one Pion only

      if(nPi0 > 1 && numuCCpi0 == true) {
          //std::cout << "Multiple PiZeros in MC: " << nPi0 << std::endl;
      }
      //else return false;
      
      
      
    
      
      // now only looking at numuCCpi0 flagged events
      
      if(numuCCpi0) {               // now already know this event has at least Pi0, one muon
          
          // Vertex
          TVector3 vtx = {-9999,-9999,-9999};
          
          for(auto mcp : all_particles) {                                                                                                       //iterate through MCParticles
              
                      if(mcp.PdgCode() == 111 && mcp.StatusCode() == 1 && foundOneMC == false) {                                                // find pizero
                          
                          int n_s = 0;
                          int n_s_photon = 0;
                          int n_m = 0;
                          int n_t = 0;
                          double muon_len = 0;
                          
                          //count tracks
                          for(auto mctrk : *mctrk_v) {
                              if(abs(mctrk.PdgCode()) == 13 && mctrk.Process() == "primary") {
                                  if(length(mctrk) > 0) {
                                      vtx.SetXYZ(mctrk.Start().X(), mctrk.Start().Y(), mctrk.Start().Z());                              // set vertex to be the start point of the primary muon
                                      n_m++;                                                                                            // primary muons
                                      muon_len = length(mctrk);
                                  }
                              }
                              else if(mctrk.PdgCode() != 2112 && mctrk.PdgCode() != 130 && mctrk.PdgCode() != 310 && mctrk.PdgCode() != 111
                                      && mctrk.PdgCode() != 311 && mctrk.PdgCode() < 10000 && mctrk.Process() == "primary" && length(mctrk) != 0 && StartInTPC(mctrk)) {
                                  n_t++;                                                                                                                  // other tracks from interaction
                              }
                          }
                          
                          //count showers
                          for(auto mcshw : *mcshow_v) {
                              if(mcshw.PdgCode() == 22 && StartInTPC(mcshw) && mcshw.MotherPdgCode() == 111) {                  //count photon showers from pizero
                                  n_s_photon++;
                              }
                              else if((mcshw.PdgCode() == 22 || abs(mcshw.PdgCode()) == 11) && StartInTPC(mcshw)) {             // count other showers
                                  TVector3 svtx;
                                  svtx.SetXYZ(mcshw.Start().X(), mcshw.Start().Y(), mcshw.Start().Z());
                                  if((svtx - vtx).Mag() < 1) {
                                      n_s++;
                                  }
                              }
                          }
                          
                          
                          // assess topology
                          //std::cout << "\t n_s " << n_s << "\t n_m " << n_m << "\t n_t " << n_t << std::endl;
                          
                          if(n_s_photon == 2 && n_m == 1 && n_t >=1 && n_s == 0) {                                                                             // selection Criteria
                              topoCut = true;
                              foundOneMC = true;
                          }
                      }
              
              // if there is already a MC pion interaction that passes:
                else if(mcp.PdgCode() == 111 && mcp.StatusCode() == 1 && foundOneMC == true) {
                    
                    int n_s = 0;
                    int n_s_photon = 0;
                    int n_m = 0;
                    int n_t = 0;
                    double muon_len = 0;
                    
                    //count tracks
                    for(auto mctrk : *mctrk_v) {
                        if(abs(mctrk.PdgCode()) == 13 && mctrk.Process() == "primary") {
                            if(length(mctrk) > 0) {
                                vtx.SetXYZ(mctrk.Start().X(), mctrk.Start().Y(), mctrk.Start().Z());                              // set vertex to be the start point of the primary muon
                                n_m++;                                                                                            // primary muons
                                muon_len = length(mctrk);
                            }
                        }
                        else if(mctrk.PdgCode() != 2112 && mctrk.PdgCode() != 130 && mctrk.PdgCode() != 310 && mctrk.PdgCode() != 111       // other tracks from interaction
                                && mctrk.PdgCode() != 311 && mctrk.PdgCode() < 10000 && mctrk.Process() == "primary" && length(mctrk) != 0 && StartInTPC(mctrk)) n_t++;
                        
                    }
                    
                    //count showers
                    for(auto mcshw : *mcshow_v) {
                        if(mcshw.PdgCode() == 22 && StartInTPC(mcshw) && mcshw.MotherPdgCode() == 111) n_s_photon++;                 //count photon showers from pizero
                            
                        else if((mcshw.PdgCode() == 22 || abs(mcshw.PdgCode()) == 11) && StartInTPC(mcshw)) {                        // count other showers
                            TVector3 svtx;
                            svtx.SetXYZ(mcshw.Start().X(), mcshw.Start().Y(), mcshw.Start().Z());
                            if((svtx - vtx).Mag() < 1) n_s++;
                        }
                    }
                    if(n_s_photon == 2 && n_m == 1 && n_t >=1 && n_s == 0) {
                        moreThanOneNuMC = true;
                    }
                    
                    
                }
              
          }
        
    } // selecting numuCCpi0 events

      if(p_topoCut && numuCCpi0 == false && topoCut == false && moreThanOneNu == false && moreThanOneNuMC == false) fakeNotCC++;
      if(p_topoCut && numuCCpi0 && topoCut == false && moreThanOneNu == false && moreThanOneNuMC == false) fakeNotTopo++;
    
      if(topoCut == true && moreThanOneNuMC == false) nTopoCut++;

      
      if(topoCut == true && p_topoCut == true && moreThanOneNuMC == false && moreThanOneNu == false) nSuccess++;
      if(topoCut == true && p_topoCut == true && moreThanOneNuMC == true && moreThanOneNu == true) nSuccess++;
      if(p_topoCut && topoCut == false && moreThanOneNu == false && moreThanOneNuMC == false) nFake++;
      
      if(moreThanOneNuMC) nMoreThanOne++;
      
      
    return true;
  }
    
//=================================================================================================//
//=================================================================================================//

  bool MCRecoCompare::finalize() {
      
      std::cout << "\n number of all events processed: \t\t" << all_events << std::endl;
      std::cout << " number of CC events \t\t\t\t" << CC_events << std::endl;
      std::cout << " number of numuCCPi0 events \t\t\t" << numuCCpi0_events << std::endl;
      std::cout << " number passing MC topoCut  \t\t\t" << nTopoCut << std::endl;
      std::cout << " number passing p_topoCut \t\t\t" << p_nTopoCut << std::endl;
      
      std::cout << " MC topoCut true and p_topoCut true \t\t" << nSuccess << std::endl;
      std::cout << " MC topoCut false and p_topoCut true (Fake) \t" << nFake << std::endl;
      std::cout << " \t of which not numuCCpi0 \t\t" << fakeNotCC << std::endl;
      std::cout << " \t of which not correct topology in MC \t" << fakeNotTopo << std::endl;
      std::cout << "\n \t more than one " << nMoreThanOne   << std::endl;
      
    
      
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
