#ifndef LARLITE_MCRECOCOMPARE_CXX
#define LARLITE_MCRECOCOMPARE_CXX

#include "MCRecoCompare.h"

namespace larlite {

  bool MCRecoCompare::initialize() {

  // FOR MC
      
      CCpi0_events      = 0;
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
      

  // FOR RECO
      
      p_CCpi0_events     = 0;
      p_all              = 0;
      p_nTopoCut         = 0;
      
   

    return true;
  }
  
//=================================================================================================//
//=================================================================================================//
    
  bool MCRecoCompare::analyze(storage_manager* storage) {
  

      std::cout << "\n ======================================================== " <<std::endl;
      std::cout << " EVENT NUMBER " << storage->event_id() << std::endl;
      
      // Reco data products
      auto pfp_v = storage->get_data<event_pfpart>("pandoraNu");
      
      std::cout << pfp_v->at(0).size() << std::endl;
      
    
      
      
      //MCTruth data products
      auto mc = storage->get_data<event_mctruth>("generator");
      auto& all_particles = mc->at(0).GetParticles();
      auto& nu = mc->at(0).GetNeutrino();
      
      int nMu       = 0;        //number of MCTruth muons
      int nPi0      = 0;        // number of MCTruth Pi0
      
      bool numuCCpi0 = false;   //flag for numuCCpi0 event
      bool topoCut   = false;   //flag for topology cut
      
      
      
      //MCTrack / MCShower data products
      
      int nMuTrk    = 0;        //number of primary muon tracks
      int nTrk      = 0;        //number of tracks other than primary muons (may include non-primary muons)
      int prim_trk  = 0;        //number of primary tracks other than muons (what else comes out of int)
      int nShow     = 0;        //number of showers
      int prim_show = 0;        //number of showers from the pi0 decay
      
      
      //MCTracks and MCShowers
      auto mctrk_v = storage->get_data<event_mctrack>("mcreco");  // this is a vector of tracks
      auto mcshow_v = storage->get_data<event_mcshower>("mcreco");// this is a vector of showers
      
      
      // Check whether event is numuCCpi0
      if(nu.CCNC() == 0 && abs(nu.Nu().PdgCode()) == 14) {
          CCpi0_events++;
          for(auto const& par: all_particles) {
              if(par.StatusCode() == 1) {                                                                           //look at stable final states only
                  if(par.PdgCode() == 111)      nPi0++;
                  if(abs(par.PdgCode()) == 13)  nMu++;
              }
          }
      } // select numu CC event
      
      
      //these are the events we want to initially select and check their tracks and showers
      if(nPi0 >= 1 && nMu == 1) {
          numuCCpi0 = true;
          all++;
      }
      
      //////////////////////////////////////////////////////////////////////////////////
      // This bit only executed when looking at a numuCCpi0 event                     //
      
      if(numuCCpi0 == true) {
          
          // Vertex
          TVector3 vtx = {-9999,-9999,-9999};
          
          for(auto mctrk : *mctrk_v){
              if(StartInTPC(mctrk) == true){                                                                        //ensure looking at tracks inside TPC volume
                  
                  
                  if(abs(mctrk.PdgCode()) == 13 && mctrk.Process() == "primary") {
                      nMuTrk++;                                                                                         // primary muon count
                      TVector3 start  = { mctrk.Start().X(), mctrk.Start().Y(), mctrk.Start().Z()};
                      TVector3 end    = { mctrk.End().X(),mctrk.End().Y(),mctrk.End().Z() };
                      muon_length = length(mctrk);                                                                      // using length function - mcsteps added up
                      vtx = start;                                                                                      // vtx: start of primary muon,from numu interaction
                      std::cout << " \n vertex at " << vtx.X() << ","<<vtx.Y()<<","<<vtx.Z()<<std::endl;
                      std::cout << " muon length " << muon_length << std::endl;
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
                      std::cout << "\n other track in event " << mctrk.Start().X()<<", "<< mctrk.Start().Y()<<", " << mctrk.Start().Z()<< std::endl;
                      std::cout << " other track length " << other_length << std::endl;
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
                      std::cout << " \n shower from " <<mcshow.Start().X()<<","<<mcshow.Start().Y()<<","<<mcshow.Start().Z()<< std::endl;
                  }
              }
          } // acting on each element of mcshow_v
          
          // count events that pass the topoCut
          if(prim_trk == 1 && nMuTrk == 1 && prim_show == 2) {
              topoCut = true;                                                   // topoCut is 1 muon track, one other track and 2 showers from vertex
              nTopoCut++;
          }
          
      } // selecting numuCCpi0 events
      
      //////////////////////////////////////////////////////////////////////////////////////////
      
      lengthRatio = other_length/muon_length;
      
      //Filling appropriate trees ==== NEEDS WORK
      if(topoCut) {
          muons->Fill();
          other->Fill();
          tree->Fill();
      }


  
    return true;
  }
    
//=================================================================================================//
//=================================================================================================//

  bool MCRecoCompare::finalize() {
      
      std::cout << "\n number of CCPi0 events \t" << CCpi0_events << std::endl;
      std::cout << " number of numuCCPi0 events \t" << all << std::endl;
      std::cout << " number passing topoCut  \t" << nTopoCut << std::endl;
      
      
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
