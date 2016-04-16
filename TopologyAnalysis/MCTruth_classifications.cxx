#ifndef LARLITE_MCTRUTH_CLASSIFICATIONS_CXX
#define LARLITE_MCTRUTH_CLASSIFICATIONS_CXX

#include "MCTruth_classifications.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  // This is run once at the beginning of the job
  bool MCTruth_classifications::initialize() {

    //make sure you initialize all your variables here:
    all_events = 0;
    CC_events = 0;
    CCpi0_events = 0; 
    
    //Add new elements to the TTree here:
    if(tree) delete tree;
    tree = new TTree("tree","tree");
    tree->Branch("thisInt",&thisInt,"thisInt/I");
    tree->Branch("thisDouble",&thisDouble,"thisDouble/D");

    return true;
  }
  
  //This is run once for every "ARTEvent", basically just 1 readout
  bool MCTruth_classifications::analyze(storage_manager* storage) {
    
    //If you want to use MCShowers and MCTracks you can use lines like this
    //
    // auto mctrk_v = storage->get_data<event_mctrack>("mcreco");
    // auto mcshow_v = storage->get_data<event_mcshower>("mcreco");
    // 
 
    auto mc = storage->get_data<event_mctruth>("generator");

    //vector of mcpart 
    auto& all_particles = mc->at(0).GetParticles();
    
    //single mcnu
    auto& nu = mc->at(0).GetNeutrino();
    
    
    // If interaction is CC and numu
    if(nu.CCNC() == 0 && abs(nu.Nu().PdgCode()) == 14){
      
      // iterate through all the particles
      // remember that these particles are 
      // at GENIE level, this means that you
      // check "how long they are" instead
      // you would need to use MCTracks for
      // that type of study
      for(auto const& par: all_particles){
	
	/// SUPER ANNOYING cout just for you :) 
	std::cout << par.PdgCode() << std::endl; 
	
      }// iterate over all particles
    }// select numu CC events
  
      
    
    
    tree->Fill();
    
    return true;
  }
  
  //This is run once at the end of the job
  bool MCTruth_classifications::finalize() {
    
    if(_fout) { _fout->cd(); tree->Write(); }
    
    return true;
  }

}
#endif
