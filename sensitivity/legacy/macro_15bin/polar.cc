const int ncosBins = 9;
const int nEBins = 10;

const double alpha = 0.65; //Asymmetry parameter of Lambda
const double GeVToMeV = 1000.;
const double M_Proton = 938.272;


const double p_min = 720.;
const double p_max = 770.;

void polar(){
  TFile* file = new TFile("../data/7201_7202_CSFlat.root");
  TTree* tree = (TTree*)file->Get("g4hyptpc");

  vector<TParticle> *BEAM = nullptr;
  vector<TParticle> *PRM = nullptr;
  vector<TParticle> *SEC = nullptr;

  double cos_theta;

  TBranch *b_BEAM = tree->GetBranch("BEAM");
  b_BEAM->SetAddress(&BEAM);

  TBranch *b_PRM = tree->GetBranch("PRM");
  b_PRM->SetAddress(&PRM);

  TBranch *b_SEC = tree->GetBranch("SEC");
  b_SEC->SetAddress(&SEC);

  tree->SetBranchAddress("cos_theta",&cos_theta);
  

  TH1D *hist_plambda = new TH1D("hist_plambda","hist_plambda",1000,-10,10);
  
  TH2D *hist_Pol = new TH2D("hist_pol","hist_pol",ncosBins,-1,1,ncosBins,-10,10);

  double d_sig_eppi[ncosBins];
  double d_sig_eL[ncosBins];

  double N_acc[nEBins][ncosBins];
  for(int i=0;i<ncosBins;i++){
    d_sig_eppi[i] = 0.;
    d_sig_eL[i] = 0.;
  }
  
  for(int i=0;i<nEBins;i++){
    for(int j=0;j<ncosBins;j++){
      N_acc[i][j] = 0.;
    }
  }
  
  for(int n=0;n<tree->GetEntries();n++){
    //for(int n=0;n<200000;n++){

    tree->GetEntry(n);

    
    const auto &p_BEAM = (*BEAM)[0];
    
    TVector3 p_Kaon(p_BEAM.Px(),p_BEAM.Py(),p_BEAM.Pz());
    double beam_mom = p_Kaon.Mag();
    TLorentzVector LVKaon(p_Kaon, TMath::Hypot(p_Kaon.Mag(), p_BEAM.GetMass()*GeVToMeV));
    TLorentzVector LVTarget(0., 0., 0. , M_Proton);
    TLorentzVector LVW = LVKaon + LVTarget;
    TVector3 beta_W = -LVW.BoostVector();
    
    TLorentzVector LVKaon_CM = LVKaon;
    TLorentzVector LVTarget_CM = LVTarget;
    LVKaon_CM.Boost(beta_W);
    LVTarget_CM.Boost(beta_W);

    //Decay
    TLorentzVector LVLambda_CM;
    TLorentzVector LVeta_CM;
    for(const auto &p_PRM : *PRM){
      if(p_PRM.GetPdgCode() == 3122){
	TVector3 p_L(p_PRM.Px(), p_PRM.Py(), p_PRM.Pz());
	LVLambda_CM.SetPxPyPzE(p_PRM.Px(), p_PRM.Py(), p_PRM.Pz(),TMath::Hypot(p_L.Mag(),p_PRM.GetMass()*GeVToMeV));
      }
      else if(p_PRM.GetPdgCode() == 221){
	TVector3 p_eta(p_PRM.Px(), p_PRM.Py(), p_PRM.Pz());
	LVeta_CM.SetPxPyPzE(p_PRM.Px(), p_PRM.Py(), p_PRM.Pz(),TMath::Hypot(p_eta.Mag(),p_PRM.GetMass()*GeVToMeV));
      }
    }
    
    LVLambda_CM.Boost(beta_W);
    LVeta_CM.Boost(beta_W);

    //Reaction plane
      
    //Unit MeV
    vector<TLorentzVector> LVDecays;
    TLorentzVector LVProton;
    TLorentzVector LVLambda;
    
    bool L_p = false;
    bool L_pi = false;
    for(const auto &p_SEC : *SEC){
      if(p_SEC.GetMother(0) == 3122){
	TVector3 p_decay(p_SEC.Px(), p_SEC.Py(), p_SEC.Pz());
	TLorentzVector LVDecay(p_decay, TMath::Hypot(p_decay.Mag(),p_SEC.GetMass()*GeVToMeV));
	LVDecays.push_back(LVDecay);
	if(p_SEC.GetPdgCode() == 2212){
	  LVProton = LVDecay;
	  L_p = true;
	}
	if(p_SEC.GetPdgCode() == -211)L_pi = true;
      }
    }
    for(size_t i=0;i<LVDecays.size();i++){
      LVLambda += LVDecays[i];
    }

    int cos_bin = std::floor((cos_theta + 1.0) / (2. / (double)ncosBins));
    int E_bin = std::floor((beam_mom - p_min) / ((p_max - p_min) / (double)nEBins));
    
  
    
    //Boost to Lambda rest frame
    if(L_p && L_pi){
      TVector3 beta_L = -LVLambda.BoostVector();
      LVProton.Boost(beta_L);

      LVKaon.Boost(beta_L);
      LVTarget.Boost(beta_L);

      TVector3 p_Kaon_L = LVKaon.Vect();
      TVector3 p_Target_L = LVTarget.Vect();

      TVector3 p_L = p_Kaon_L.Cross(p_Target_L);

      //Get Proton Direction Vector from Lambda rest frame
      TVector3 p_Proton_L = LVProton.Vect();
      double cosThetaP = p_Proton_L.Unit().Dot(p_L.Unit());

      d_sig_eppi[cos_bin] += cosThetaP;
    }
    d_sig_eL[cos_bin] +=1.;

    hist_plambda->Fill(LVLambda.P());


    //Acceptance per Energy / Cos theta (eta)
    N_acc[E_bin][cos_bin] += 1.;
    
  }
  cout<<"1"<<endl;

  TGraph *pol_L = new TGraph();
  TGraph *graph_acc = new TGraph();
  for(int i=0;i<ncosBins;i++){
    double step = 2. / (double)ncosBins;
    double cos_eta = -1 + step/2. + step*i;
    cout<<(3.0 / alpha) * d_sig_eppi[i] / d_sig_eL[i]<<endl;
    pol_L->SetPoint(pol_L->GetN(),cos_eta,(3.0 / alpha) * d_sig_eppi[i] / d_sig_eL[i]);
    graph_acc->SetPoint(graph_acc->GetN(),cos_eta,N_acc[5][i]);
  }

  cout<<"2"<<endl;
  



  TCanvas *c1 = new TCanvas("c1","c1");
  //hist_plambda->Draw();
  //hist_plambda_prm->SetLineColor(kRed);
  //hist_plambda_prm->Draw("same");
  //hist_Pol->Draw();
  pol_L->SetMarkerStyle(20);
  pol_L->GetYaxis()->SetRangeUser(-1.0, 1.0);
  pol_L->Draw("AP");


  cout<<"3"<<endl;


  /*
  TCanvas *c2 = new TCanvas("c2","c2");
  graph_acc->SetMarkerStyle(20);
  graph_acc->Draw("AP");
  */


  //gApplication->Terminate();
}
