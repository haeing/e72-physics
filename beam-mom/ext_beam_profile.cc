#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

const double target_pos = -143 - 150;
const double start_pos = -900.; //Geant4 coordinate is different!! This position is E72 experimental coordinate 

void ext_beam_profile(){

  int mom = 735;
  int runnumber = 2447;
  
  TString dir = "/gpfs/group/had/sks/Users/haein/data/JPARC2025Nov_root/blc";
  TFile *file_track = new TFile(Form("%s/run0%d_K18Tracking.root",dir.Data(),runnumber)); //Beam
  TTree *tree_track = (TTree*)file_track->Get("k18");
  

  TFile *file_beam = new TFile(Form("/hsm/had/sks/E72/JPARC2025Nov/beam_simul/beam_profile_run0%d_-%d.root",runnumber,mom), "RECREATE");
  TTree * tree_beam = new TTree("tr","Beam Profile");


  std::vector<double> *px = nullptr;
  std::vector<double> *py = nullptr;
  std::vector<double> *pz = nullptr;
  std::vector<double> *xout = nullptr;
  std::vector<double> *yout = nullptr;
  std::vector<double> *uout = nullptr;
  std::vector<double> *vout = nullptr;
  
  
  double pointInx,pointIny,pointInz,pInx,pIny,pInz;

  tree_track->SetBranchAddress("px",&px);
  tree_track->SetBranchAddress("py",&py);
  tree_track->SetBranchAddress("pz",&pz);
  tree_track->SetBranchAddress("xout",&xout);
  tree_track->SetBranchAddress("yout",&yout);
  tree_track->SetBranchAddress("uout",&uout);
  tree_track->SetBranchAddress("vout",&vout);


  tree_beam->Branch("pointInx",&pointInx,"pointInx/D");
  tree_beam->Branch("pointIny",&pointIny,"pointIny/D");
  tree_beam->Branch("pointInz",&pointInz,"pointInz/D");
  tree_beam->Branch("pInx",&pInx,"pInx/D");
  tree_beam->Branch("pIny",&pIny,"pIny/D");
  tree_beam->Branch("pInz",&pInz,"pInz/D");
  


  for(int i=0;i<tree_track->GetEntries();i++){
    tree_track->GetEntry(i);

    const auto ntrack = std::min({px->size(), py->size(), pz->size(),
                                  xout->size(), yout->size(),
                                  uout->size(), vout->size()});

    for(size_t itrack=0; itrack<ntrack; ++itrack){
      if(pz->at(itrack)<0)continue;

      pointInx = xout->at(itrack)-start_pos*uout->at(itrack);
      pointIny = yout->at(itrack)-start_pos*vout->at(itrack);
      pointInz = start_pos + 150.071;

      pInx = px->at(itrack);
      pIny = py->at(itrack);
      pInz = pz->at(itrack);

      tree_beam->Fill();
    }
  }
  
  std::cout<<"write close"<<std::endl;

  tree_beam->Write();
  file_beam->Close();
 
}
