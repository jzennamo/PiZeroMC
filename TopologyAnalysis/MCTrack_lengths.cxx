#ifndef LARLITE_MCTRACK_LENGTHS_CXX
#define LARLITE_MCTRACK_LENGTHS_CXX

#include "MCTrack_lengths.h"
#include "DataFormat/mctruth.h"
#include <TLorentzVector.h>

namespace larlite {

  bool MCTrack_lengths::initialize() {

    //
    // This function is called in the beggining of event loop
    // Do all variable initialization you wish to do here.
    // If you have a histogram to fill in the event loop, for example,
    // here is a good place to create one on the heap (i.e. "new TH1D"). 
    //
      
      CCpi0_events      = 0;
      all               = 0; //these are all CCpi0 events with one muon, one pi0 (numuCCPi0)
      nTopoCut          = 0; //number of events passing topoCut
      
      if(tree) delete tree;
      tree = new TTree("tree","tree");
      tree->Branch("thisInt",&thisInt,"thisInt/I");
      tree->Branch("no_muons",&no_muons,"no_muons/I");
      tree->Branch("no_tracks",&no_tracks,"no_tracks/I");                   // number of tracks per event other than muons
      tree->Branch("no_prim_tracks",&no_prim_tracks,"no_prim_tracks/I");
      
      muons = new TTree("muons","muons");
      muons->Branch("muon_length",&muon_length,"muon_length/D");
      
      other = new TTree("other","other");
      other->Branch("other_length",&other_length,"other_length/D");
      other->Branch("other_pdg",&other_pdg,"other_pdg/D");

    return true;
  }
  
  //This is run once for every "ARTEvent"
  bool MCTrack_lengths::analyze(storage_manager* storage) {
  
      
    //
    // Do your event-by-event analysis here. This function is called for 
    // each event in the loop. You have "storage" pointers which contains
    // event-wise data. To see what is available, check the "Manual.pdf":
    // Or you can refer to Base/DataFormatConstants.hh for available data type
    // enum values. Here is one example of getting PMT waveform collection.
    //
      
      //MCTruth data products
      auto mc = storage->get_data<event_mctruth>("generator");
      auto& all_particles = mc->at(0).GetParticles();
      auto& nu = mc->at(0).GetNeutrino();
    
      int nMu       = 0;        //number of MCTruth muons
      int nPi0      = 0;        // number of MCTruth Pi0
      
      bool numuCCpi0 = false;   //flag for numuCCpi0 event
      bool topoCut   = false;   //flag for topology cut
      
      
      //MCTrack / MCShower data products
      
      int nMuTrk    = 0;        //number of muon tracks
      int nTrk      = 0;        //number of tracks other than muons
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
            
          if(par.StatusCode() == 1) {
              
              if(par.PdgCode() == 111)      nPi0++;
              if(abs(par.PdgCode()) == 13)  nMu++;
              
          } //StatusCode() == 1 (stable final state particels)
        } // iterate over all MCTruth particles
      } // select numu CC event
      
      
      //these are the events we want to initially select and check their tracks and showers
      if(nPi0 >= 1 && nMu == 1) {
          numuCCpi0 = true;
          all++;
      }
      
      //////////////////////////////////////////////////////////////////////////////////
      // This bit only executed when looking at a numuCCpi0 event                     //
      
      if(numuCCpi0 == true) {     // these are the events that contribute to all count
          //std::cout << "\n\n" <<std::endl;
          
        // Vertex
        TVector3 vtx;
          
        for(auto mctrk : *mctrk_v){
            
          if(abs(mctrk.PdgCode()) == 13) {
              
            nMuTrk++;
            double mu_len = 0; TVector3 start; TVector3 end;
            start.SetXYZ(mctrk.Start().X(),mctrk.Start().Y(),mctrk.Start().Z());
            end.SetXYZ(mctrk.End().X(),mctrk.End().Y(),mctrk.End().Z());
            mu_len = (end - start).Mag();
            muon_length = mu_len;
            muons->Fill();
              
            // Define vertex as start point of "primary" muon, i.e. the one coming from numu interaction
              if(mctrk.Process() == "primary") {
                  vtx = start ;
              }
              else vtx = {-9999,-9999,-9999};
              //std::cout << "process " << mctrk.Process() << std::endl;
              //std::cout << "vtx \t"      << vtx.X() << std::endl;
          }
          
          
          else if(mctrk.PdgCode() != 2112 && mctrk.PdgCode() != 111 && mctrk.PdgCode() < 10000) {            // cutting out neutral particles and fragments
            nTrk++;
            double other_len = 0; TVector3 start; TVector3 end;
            start.SetXYZ(mctrk.Start().X(),mctrk.Start().Y(),mctrk.Start().Z());
            end.SetXYZ(mctrk.End().X(),mctrk.End().Y(),mctrk.End().Z());
            other_len = (end - start).Mag();
            other_length = other_len;
            d_vtx = (vtx - start).Mag(); // distance to vertex
            other_pdg = mctrk.PdgCode(); //std::cout << other_pdg << std::endl;
            other_process = mctrk.Process();
            if(other_process == "primary") prim_trk++;
            other->Fill();
          }
            
            //std::cout << "\n distance to vertex " << d_vtx << std::endl;
            //std::cout << " other pdg " << other_pdg << std::endl;
            //std::cout << " other process " << other_process << std::endl;
        
        } //acting on each element of mctrk_v
          
          for(auto mcshow : *mcshow_v) {
              nShow++;
              if(mcshow.MotherPdgCode() == 111 && mcshow.MotherTrackID() == mcshow.AncestorTrackID()) prim_show++;
              //std::cout << "shower process " << mcshow.Process() << std::endl;
              //std::cout << "shower mother " << mcshow.MotherPdgCode() << std::endl;
              //std::cout << "mother track ID " << mcshow.MotherTrackID() <<std::endl;
              //std::cout << "ancestor track ID " << mcshow.AncestorTrackID() << std::endl;
          }
          
          if(prim_trk == 1 && nMuTrk == 1 && prim_show == 2) {
              topoCut = true;   // topoCut is 1 muon track, one other track and 2 showers from vertex
              nTopoCut++;
          }
          
          //std::cout << "topoCut " << topoCut << std::endl;
          //std::cout << "nShow " << nShow << std::endl;
    } // selecting numuCCpi0 events
      
    //////////////////////////////////////////////////////////////////////////////////////////
     
    
    //if(muon_length > 10000) std::cout << "muon_length "<< muon_length << std::endl;
      
      
    
    no_muons = nMuTrk;
    if(numuCCpi0 == false) no_muons = -3; // -3 if no interaction present
    no_tracks = nTrk;
    if(numuCCpi0 == false) no_tracks = -3; // -3 if no interaction present
    no_prim_tracks = prim_trk;
    if(numuCCpi0 == false) no_prim_tracks = -3; // use -3 as flag for no interaction present
    if(numuCCpi0 == true) tree->Fill();
      
    return true;
  }

  //This is run once at the end of the job
  bool MCTrack_lengths::finalize() {

      std::cout << "\n number of CCPi0 events \t" << CCpi0_events << std::endl;
      std::cout << " number of all           \t" << all << std::endl;
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


// Reference larlite manual:
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
