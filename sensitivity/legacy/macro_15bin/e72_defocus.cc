#include <TParticle.h>
#include <vector>
#include "../include/TPCPadHelper.hh"
//#inlcude "/Users/ihaein/Work/E72/HypTPC/TPCPadHelpher.hh"

void e72_defocus(){
  //TH2Poly *TPC2DPoly = new TH2Poly("TPC2DPoly", "TPC2DPoly", MinZ, MaxZ, MinX, MaxX);
  TH2Poly *TPC2DPoly = new TH2Poly("TPC2DPoly", "TPC2DPoly", -300, 300, -300, 300);

  TH2D *hist_yz = new TH2D("hist_yz","hist_yz",300,-300,300,300,-300,300);

  Double_t X[5];
  Double_t Y[5];
  for (Int_t l=0; l<NumOfLayersTPC; ++l) {
    Double_t pLength = padParameter[l][5];
    Double_t st      = (180.-(360./padParameter[l][3]) *
			padParameter[l][1]/2.);
    Double_t sTheta  = (-1+st/180.)*TMath::Pi();
    Double_t dTheta  = (360./padParameter[l][3])/180.*TMath::Pi();
    Double_t cRad    = padParameter[l][2];
    Int_t    nPad    = padParameter[l][1];
    for (Int_t j=0; j<nPad; ++j) {
      Bool_t dead = false;
      for(int i=0;i<sizeof(padOnFrame)/sizeof(*padOnFrame);i++){
	if(GetPadId(l,j) == padOnFrame[i]){
	  dead=true;
	}
	
      }

      if(!dead){
	X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
	X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
	X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
	X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
	X[0] = X[4];
	Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
	Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
	Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
	Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
	Y[0] = Y[4];
	for (Int_t k=0; k<5; ++k) X[k] += ZTarget;
      }
      TPC2DPoly->AddBin(5, X, Y);

    }
  }

  TFile *fFile_g4 = new TFile("data/e72_defocus_HS_off.root");
  TTree *tree = (TTree*)fFile_g4->Get("g4hyptpc");

  vector<TParticle> *TPC = nullptr;
  TBranch *b_TPC = tree->GetBranch("TPC");
  b_TPC->SetAddress(&TPC);

  int nhittpc;
  double xtpc[5000];
  double ytpc[5000];
  double ztpc[5000];

  tree->SetBranchAddress("nhittpc",&nhittpc);
  tree->SetBranchAddress("xtpc",xtpc);
  tree->SetBranchAddress("ytpc",ytpc);
  tree->SetBranchAddress("ztpc",ztpc);
  
  
  for(int n=0;n<tree->GetEntries();n++){
    tree->GetEvent(n);

    for(int i=0;i<nhittpc;i++){
      TPC2DPoly->Fill(ztpc[i],xtpc[i]);

      hist_yz->Fill(ztpc[i],ytpc[i]);
    }
    /*
    for(const auto &p_TPC : *TPC){
      TPC2DPoly->Fill(p_TPC.Pz(),p_TPC.Px());
    }
    */
    
    

  }

  TCanvas *c1 = new TCanvas("c1","c1");
  gStyle->SetOptStat(0);
  TPC2DPoly->SetTitle(";X [mm];Z [mm]");
  TPC2DPoly->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","c2");
  hist_yz->Draw("colz");
 
}
