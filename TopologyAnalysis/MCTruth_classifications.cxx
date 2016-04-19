#ifndef LARLITE_MCTRUTH_CLASSIFICATIONS_CXX
#define LARLITE_MCTRUTH_CLASSIFICATIONS_CXX

#include "MCTruth_classifications.h"
#include "DataFormat/mctruth.h"

namespace larlite {

  // This is run once at the beginning of the job
  bool MCTruth_classifications::initialize() {

    //make sure you initialize all your variables here:
    all_events          = 0;
    CC_events           = 0;
    CCpi0_events        = 0;
    all                 = 0;
    nothing             = 0;
    nothing_n           = 0;
    single_proton       = 0;
    single_proton_n     = 0;
    single_chPi         = 0;
    single_chPi_n       = 0;
    single_pi0          = 0;
    single_pi0_n        = 0;
    proton_pi0          = 0;
    proton_pi0_n        = 0;
    proton_chPi         = 0;
    proton_chPi_n       = 0;
    proton_more2        = 0;
    proton_more2_n      = 0;
    chPi_more2          = 0;
    chPi_more2_n        = 0;
    chPi_more2_pros     = 0;
    chPi_more2_pros_n   = 0;
    Pi0_more2pro        = 0;
    Pi0_more2pro_n      = 0;
    chPi_more2pro       = 0;
    chPi_more2pro_n     = 0;
    more_chPi_pro_pi0   = 0;
    more_chPi_pro_pi0_n = 0;
    least1Kaon          = 0;
    morePi0_pros        = 0;
    morePi0_pros_n      = 0;
    nuclear             = 0;
    
    
    
    //Add new elements to the TTree here:
    if(tree) delete tree;
    tree = new TTree("tree","tree");
    tree->Branch("thisInt",&thisInt,"thisInt/I");
    tree->Branch("thisDouble",&thisDouble,"thisDouble/D");
    tree->Branch("fpdg",&fpdg,"fpdg/I");
    tree->Branch("nparticles",&nparticles,"nparticles/I");                  // number of final state particles in the interaction
    tree->Branch("fProtons",&fProtons,"fProtons/I");                  // number of final state particles in the interaction that are protons
      
    return true;
  }
  
  //This is run once for every "ARTEvent", basically just 1 readout
  bool MCTruth_classifications::analyze(storage_manager* storage) {
     
    all_events++;
    bool containsPiZero = false;
    bool containsMuPiZero = false;
    
    //If you want to use MCShowers and MCTracks you can use lines like this
    //
    // auto mctrk_v = storage->get_data<event_mctrack>("mcreco");
    // auto mcshow_v = storage->get_data<event_mcshower>("mcreco");
    // 
 
    auto mc = storage->get_data<event_mctruth>("generator");
      
    //std::cout<<"\n number of particles in event: " << mc->at(0).NParticles() << std::endl;
      
    //vector of mcpart
    auto& all_particles = mc->at(0).GetParticles();
    
    //single mcnu
    auto& nu = mc->at(0).GetNeutrino();
    
    int all_final_stable    = 0;
    int nMu                 = 0;
    int nProtons            = 0;
    int nPi0                = 0;
    int nPiCh               = 0;
    int nNeutrons           = 0;
    int nKaons              = 0;
    int nFrag               = 0;

    
    // If interaction is CC and numu
    if(nu.CCNC() == 0 && abs(nu.Nu().PdgCode()) == 14){
      
      CCpi0_events++;
      containsPiZero = false;
        
      //std::cout << "\t event contains CC interaction: YES \n \t interaction type: " <<nu.InteractionType() << std::endl;
      
      // iterate through all the particles
      // remember that these particles are 
      // at GENIE level, this means that you can't
      // check "how long they are" instead
      // you would need to use MCTracks for
      // that type of study
      
      for(auto const& par: all_particles){
          
        if(par.StatusCode() == 1) {
          
	      /// SUPER ANNOYING cout just for you :)
	      //std::cout << par.PdgCode() << std::endl;
            
          all_final_stable++;
          
          //if (par.PdgCode() == 111) {
          //  containsPiZero = true;
          //}
         
          if(par.PdgCode() == 2212)                         nProtons++;
          if(par.PdgCode() == 111)                          nPi0++;
          if(par.PdgCode() == 13 )                          nMu++;
          if(par.PdgCode() == 211 || par.PdgCode() == -211) nPiCh++;
          if(par.PdgCode() == 2112)                         nNeutrons++;
          if(par.PdgCode() == 321 || par.PdgCode() == 130)  nKaons++;
          if(par.PdgCode() > 100000)                        nFrag++;
    
        
        } //StatusCode() == 1 (stable final state)
        }// iterate over all particles
      }// select numu CC events
      
    if(nMu == 1 && nPi0 >= 1){ containsMuPiZero = true;                                                                                 all++;};
    if(nMu == 1 && nPi0 == 1 && all_final_stable == 2)                                                                                  nothing ++;
    else if(nMu == 1 && nPi0 == 1 && nNeutrons > 0 && all_final_stable == nNeutrons+2)                                                  nothing_n++;
    else if(nMu == 1 && nPi0 == 1 && nProtons == 1 && all_final_stable == 3)                                                            single_proton++;
    else if(nMu == 1 && nPi0 == 1 && nProtons == 1 && nNeutrons > 0 && all_final_stable == nNeutrons+3)                                 single_proton_n++;
    else if(nMu == 1 && nPi0 == 1 && nPiCh == 1 && all_final_stable == 3)                                                               single_chPi++;
    else if(nMu == 1 && nPi0 == 1 && nPiCh == 1 && nNeutrons > 0 && all_final_stable == 3+nNeutrons)                                    single_chPi_n++;
    else if(nMu == 1 && nPi0 == 2 && all_final_stable == 3)                                                                             single_pi0++;
    else if(nMu == 1 && nPi0 == 2 && nNeutrons > 0 && all_final_stable == 3+nNeutrons)                                                  single_pi0_n++;
    else if(nMu == 1 && nPi0 == 2 && nProtons == 1 && all_final_stable == 4)                                                            proton_pi0++;
    else if(nMu == 1 && nPi0 == 2 && nProtons == 1 && nNeutrons > 0 && all_final_stable == 4+nNeutrons)                                 proton_pi0_n++;
    else if(nMu == 1 && nPi0 == 1 && nProtons == 1 && nPiCh >= 1 && all_final_stable == 3+nPiCh)                                        proton_chPi++;
    else if(nMu == 1 && nPi0 == 1 && nProtons == 1 && nPiCh >= 1 && nNeutrons > 0 && all_final_stable == 3+nPiCh+nNeutrons)             proton_chPi_n++;
    else if(nMu == 1 && nPi0 == 1 && nProtons >=2 && nPiCh == 1 && all_final_stable == 3+nProtons)                                      chPi_more2pro++;
    else if(nMu == 1 && nPi0 == 1 && nProtons >=2 && nPiCh == 1 && nNeutrons > 0 && all_final_stable == 3+nProtons+nNeutrons)           chPi_more2pro_n++;
    else if(nMu == 1 && nPi0 == 2 && nProtons >=2 && all_final_stable == 3+nProtons)                                                    Pi0_more2pro++;
    else if(nMu == 1 && nPi0 == 2 && nProtons >=2 && nNeutrons > 0 && all_final_stable == 3+nProtons+nNeutrons)                         Pi0_more2pro_n++;
    else if(nMu == 1 && nPi0 >= 3 && nProtons >= 1 && all_final_stable == 1+nPi0+nProtons)                                              morePi0_pros++;
    else if(nMu == 1 && nPi0 >= 3 && nProtons >= 1 && nNeutrons > 0 && all_final_stable == 1+nPi0+nProtons+nNeutrons)                   morePi0_pros_n++;
    else if(nMu == 1 && nPi0 == 1 && nProtons >= 2 && all_final_stable == (nProtons+2))                                                 proton_more2++;
    else if(nMu == 1 && nPi0 == 1 && nProtons >= 2 && nNeutrons > 0 && all_final_stable == nProtons+2+nNeutrons)                        proton_more2_n++;
    else if(nMu == 1 && nPi0 == 1 && nPiCh >= 2 && all_final_stable == (nPiCh+2))                                                       chPi_more2++;
    else if(nMu == 1 && nPi0 == 1 && nPiCh >= 2 && nNeutrons > 0 && all_final_stable == nPiCh+2+nNeutrons)                              chPi_more2_n++;
    else if(nMu == 1 && nPi0 == 1 && nPiCh >= 2 && nProtons >=1 && all_final_stable == 2+nPiCh+nProtons)                                chPi_more2_pros++;
    else if(nMu == 1 && nPi0 == 1 && nPiCh >= 2 && nProtons >=1 && nNeutrons > 0 && all_final_stable == 2+nPiCh+nProtons+nNeutrons)     chPi_more2_pros_n++;
    else if(nMu == 1 && nPi0 >= 2 && nProtons >= 1 && nPiCh >=1 && all_final_stable ==nPi0+nProtons+nPiCh+1)                            more_chPi_pro_pi0++;
    else if(nMu == 1 && nPi0 >= 2 && nProtons >= 1 && nPiCh >=1 && nNeutrons > 0 && all_final_stable ==nPi0+nProtons+nPiCh+1+nNeutrons) more_chPi_pro_pi0_n++;
    else if(nMu == 1 && nPi0 >= 1 && nFrag > 0)                                                                                         nuclear++;
    else if(nMu >= 1 && nPi0 >= 1 && nKaons >=1)                                                                                        least1Kaon++;
      
    else if(nMu >=1 && nPi0 >=1) std::cout << "nMu " << nMu << "\t nPi0 " << nPi0 << "\t nPiCh " << nPiCh<< "\t nProtons " << nProtons<< "\t nNeutrons " << nNeutrons << "\t nFrag " << nFrag<< "\t all" << all_final_stable << std::endl;
    
    nparticles = all_final_stable;
    fProtons = nProtons;
      
    if(containsMuPiZero) {
      tree->Fill();
    }
      
    return true;
  }
    
  
  //This is run once at the end of the job
  bool MCTruth_classifications::finalize() {
      
      std::cout << "\n\n ============================================================================= \n ";
      std::cout << " \t number of all events processed: \t\t "              << all_events       << std::endl;
      std::cout << " \t number of CCpi0 events found: \t\t\t "              << CCpi0_events     << std::endl;
      
      std::cout << "\n \t final state topolgy: muon + Pi0 + X \t\t\t\t "        << all                  << "\t\t"<< 100 << " %"<<std::endl;
      
      std::cout << " \t\t X is nothing                  \t\t\t\t "              << nothing              << "\t\t"<< 100*nothing/all << " %"<<std::endl;
      std::cout << " \t\t X is nothing and > 0 neutrons \t\t\t\t "              << nothing_n            << "\t\t"<< 100*nothing_n/all<< " %"<<std::endl;
      std::cout << " \t\t X is proton                   \t\t\t\t "              << single_proton        << "\t\t"<< 100*single_proton/all << " %"<<std::endl;
      std::cout << " \t\t X is proton and > 0 neutrons  \t\t\t\t "              << single_proton_n      << "\t\t"<< 100*single_proton_n/all<<" %"<<std::endl;
      std::cout << " \t\t X is chPi (charged Pion)      \t\t\t\t "              << single_chPi          << "\t\t"<< 100*single_chPi/all << " %"<<std::endl;
      std::cout << " \t\t X is chPi and > 0 neutrons    \t\t\t\t "              << single_chPi_n        << "\t\t"<< 100*single_chPi_n/all << " %"<<std::endl;
      std::cout << " \t\t X is Pi0                      \t\t\t\t "              << single_pi0           << "\t\t"<< 100*single_pi0/all << " %"<<std::endl;
      std::cout << " \t\t X is Pi0 and > 0 neutrons     \t\t\t\t "              << single_pi0_n         << "\t\t"<< 100*single_pi0_n/all << " %"<<std::endl;
      std::cout << " \t\t X is Pi0 + proton             \t\t\t\t "              << proton_pi0           << "\t\t"<< 100*proton_pi0/all << " %"<<std::endl;
      std::cout << " \t\t X is Pi0 + proton + > 0 neutrons \t\t\t\t "           << proton_pi0_n         << "\t\t"<< 100*proton_pi0_n/all << " %"<<std::endl;
      std::cout << " \t\t X is chPi + proton            \t\t\t\t "              << proton_chPi          << "\t\t"<< 100*proton_chPi/all << " %"<<std::endl;
      std::cout << " \t\t X is chPi + proton + > 0 neutrons \t\t\t "            << proton_chPi_n        << "\t\t"<< 100*proton_chPi_n/all << " %"<<std::endl;
      std::cout << " \t\t X is chPi + > 1 protons       \t\t\t\t "              << chPi_more2pro        << "\t\t"<< 100*chPi_more2pro/all << " %"<<std::endl;
      std::cout << " \t\t X is chPi + > 1 protons + > 0 neutrons \t\t\t "       << chPi_more2pro_n      << "\t\t"<< 100*chPi_more2pro_n/all << " %"<<std::endl;
      std::cout << " \t\t X is Pi0 + > 1 proton         \t\t\t\t "              << Pi0_more2pro         << "\t\t"<< 100*Pi0_more2pro/all << " %"<<std::endl;
      std::cout << " \t\t X is Pi0 + > 1 proton + > 0 neutrons \t\t\t "         << Pi0_more2pro_n       << "\t\t"<< 100*Pi0_more2pro_n/all << " %"<<std::endl;
      std::cout << " \t\t X is > 1 Pi0 + > 0 proton     \t\t\t\t "              << morePi0_pros         << "\t\t"<< 100*morePi0_pros/all << " %"<<std::endl;
      std::cout << " \t\t X is > 1 Pi0 + > 0 proton + > 0 neutrons \t\t\t "     << morePi0_pros_n       << "\t\t"<< 100*morePi0_pros_n/all << " %"<<std::endl;
      std::cout << " \t\t X is > 1 protons              \t\t\t\t "              << proton_more2         << "\t\t"<< 100*proton_more2/all << " %"<<std::endl;
      std::cout << " \t\t X is > 1 protons + > 0 neutrons \t\t\t\t "            << proton_more2_n       << "\t\t"<< 100*proton_more2_n/all << " %"<<std::endl;
      std::cout << " \t\t X is > 1 chPi                 \t\t\t\t "              << chPi_more2           << "\t\t"<< 100*chPi_more2/all << " %"<<std::endl;
      std::cout << " \t\t X is > 1 chPi + > 0 neutrons  \t\t\t\t "              << chPi_more2_n         << "\t\t"<< 100*chPi_more2_n/all << " %"<<std::endl;
      std::cout << " \t\t X is > 1 chPi + > 0 proton    \t\t\t\t "              << chPi_more2_pros      << "\t\t"<< 100*chPi_more2_pros/all << " %"<<std::endl;
      std::cout << " \t\t X is > 1 chPi + > 0 proton + > 0 neutrons \t\t "      << chPi_more2_pros_n    << "\t\t"<< 100*chPi_more2_pros_n/all << " %"<<std::endl;
      std::cout << " \t\t X is > 0 chPi + > 0 proton + > 0 Pi0  \t\t\t "        << more_chPi_pro_pi0    << "\t\t"<< 100*more_chPi_pro_pi0/all << " %"<<std::endl;
      std::cout << " \t\t X is > 0 chPi + > 0 proton + > 0 Pi0 + >0 neutrons \t "<< more_chPi_pro_pi0_n << "\t\t"<< 100*more_chPi_pro_pi0_n/all << " %"<<std::endl;
      std::cout << " \t\t X contains > 0 nuclear fragments \t\t\t\t "           << nuclear              << "\t\t"<< 100*nuclear/all << " %" << std::endl;
      std::cout << " \t\t X contains at least 1 kaon \t\t\t\t "                 << least1Kaon           << "\t\t"<< 100*least1Kaon/all<< " %"<<std::endl;
      
      std:: cout<< " \t\t                               \t\t\t\t"<< nothing + nothing_n + single_proton + single_proton_n + single_chPi +single_chPi_n + single_pi0 + single_pi0_n + proton_chPi + proton_chPi_n + chPi_more2pro + chPi_more2pro_n + proton_pi0 + proton_pi0_n + Pi0_more2pro + Pi0_more2pro_n + morePi0_pros + morePi0_pros_n + proton_more2 + proton_more2_n + chPi_more2 + chPi_more2_n + chPi_more2_pros + chPi_more2_pros_n + more_chPi_pro_pi0 + more_chPi_pro_pi0_n + least1Kaon + nuclear<< std::endl;
      std::cout <<" ============================================================================= \n\n\n ";
    
    if(_fout) { _fout->cd(); tree->Write(); }
    
    return true;
  }

}
#endif
