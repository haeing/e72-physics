static const double mK = 493.677;
static const double mp = 938.27208816;
static const double mn = 939.5654205;
static const double mpi = 139.57039;
static const double mpi0 = 134.9768;
static const double me = 0.51099895;

static const double meta =547.862;
static const double mLambda = 1115.683;

static const double spill_length = 4.24; // [s]
static const double intensity = 38400.;

static const double BAC_thre = 1/1.115;
static const double KVC_thre = 1/1.45;

static const double rho = 0.07085; //g/cm3
static const double NA = 6.022 *1e23; //mol-1
static const double W = 1.008; //atomic mass of LH2 target
static const double Ltarget = 7.139; //cm targe thickness-> approximately assume

static const double Ntarget = rho*NA/W*Ltarget;
double GetBeta(int pdg, double mom){
  double mass = -9999;
  if(TMath::Abs(pdg) == 321)mass = mK;
  else if(TMath::Abs(pdg) == 2212)mass = mp;
  else if(TMath::Abs(pdg) == 2112)mass = mn;
  else if(TMath::Abs(pdg) == 211)mass = mpi;
  else if(TMath::Abs(pdg) == 111)mass = mpi0;
  else if(TMath::Abs(pdg) == 11)mass = me;

  if(mass == -9999)return -9999;
  else{ return mom/TMath::Sqrt(mass*mass + mom*mom);}
  
}


void trigger_rate(){
  TFile *file1 = new TFile("../data/e72_beam.root");
  TTree *tree1 = (TTree*)file1->Get("g4hyptpc");

  TFile *file2 = new TFile("../data/7201_7207_CSOn.root");
  TTree *tree2 = (TTree*)file2->Get("g4hyptpc");

  vector<TParticle> *HTOF1 = nullptr;
  vector<TParticle> *BAC1 = nullptr;
  vector<TParticle> *KVC1 = nullptr;
  
  vector<TParticle> *HTOF2 = nullptr;
  vector<TParticle> *BAC2 = nullptr;
  vector<TParticle> *KVC2 = nullptr;
  
  TBranch *b_HTOF1 = tree1->GetBranch("HTOF");
  b_HTOF1->SetAddress(&HTOF1);

  TBranch *b_BAC1 = tree1->GetBranch("BAC");
  b_BAC1->SetAddress(&BAC1);
  
  TBranch *b_KVC1 = tree1->GetBranch("KVC");
  b_KVC1->SetAddress(&KVC1);
 
  TBranch *b_HTOF2 = tree2->GetBranch("HTOF");
  b_HTOF2->SetAddress(&HTOF2);

  TBranch *b_BAC2 = tree2->GetBranch("BAC");
  b_BAC2->SetAddress(&BAC2);
  
  TBranch *b_KVC2 = tree2->GetBranch("KVC");
  b_KVC2->SetAddress(&KVC2);

  bool mul1 = false;
  bool forward1 = false;
  bool bac1 = false;
  bool kvc1 = false;

  bool mul2 = false;
  bool forward2 = false;
  bool bac2 = false;
  bool kvc2 = false;

  int Kbeam_count = 0;
  int kp_count = 0;
  for(int n=0;n<tree2->GetEntries();n++){
    mul1 = false;
    forward1 = false;
    bac1 = false;
    kvc1 = false;

    mul2 = false;
    forward2 = false;
    bac2 = false;
    kvc2 = false;
    
    tree1->GetEntry(n);
    tree2->GetEntry(n);

    if(HTOF1->size() >=2)mul1 = true;
    for(const auto &p_HTOF1 : *HTOF1){
      if(p_HTOF1.GetMother(1)>=17 && p_HTOF1.GetMother(1)<=21 && p_HTOF1.GetWeight() > 3.0){
	forward1 = true;
      }
      
    }

    for(const auto&p_BAC1 : *BAC1){
      double BAC_pdg = p_BAC1.GetPdgCode();
      double BAC_mom = TMath::Sqrt(p_BAC1.Px()*p_BAC1.Px()+p_BAC1.Py()*p_BAC1.Py()+p_BAC1.Pz()*p_BAC1.Pz());
      double BAC_beta = GetBeta(BAC_pdg,BAC_mom);
      
      if(BAC_beta > BAC_thre){
	
	bac1 = true;
      }
    }

    if(mul1 && forward1 && !kvc1){
      Kbeam_count++;
    }
    for(const auto&p_KVC1 : *KVC1){
      double KVC_pdg = p_KVC1.GetPdgCode();
      double KVC_mom = TMath::Sqrt(p_KVC1.Px()*p_KVC1.Px()+p_KVC1.Py()*p_KVC1.Py()+p_KVC1.Pz()*p_KVC1.Pz());
      double KVC_beta = GetBeta(KVC_pdg,KVC_mom);
      
      if(KVC_beta > KVC_thre){
	
	kvc1 = true;
      }
    }

    if(mul1 || forward1){
      if(!bac1 && !kvc1){
	Kbeam_count++;
      }
    }

    //Kp->Kp
    if(HTOF2->size() >=2)mul2 = true;
    for(const auto &p_HTOF2 : *HTOF2){
      if(p_HTOF2.GetMother(1)>=27 && p_HTOF2.GetMother(1)<=22 && p_HTOF2.GetWeight() > 3.0){
	forward2 = true;
      }
      
    }

    for(const auto&p_BAC2 : *BAC2){
      double BAC_pdg = p_BAC2.GetPdgCode();
      double BAC_mom = TMath::Sqrt(p_BAC2.Px()*p_BAC2.Px()+p_BAC2.Py()*p_BAC2.Py()+p_BAC2.Pz()*p_BAC2.Pz());
      double BAC_beta = GetBeta(BAC_pdg,BAC_mom);
      
      if(BAC_beta > BAC_thre){
	
	bac2 = true;
      }
    }
    
    
    for(const auto&p_KVC2 : *KVC2){
      double KVC_pdg = p_KVC2.GetPdgCode();
      double KVC_mom = TMath::Sqrt(p_KVC2.Px()*p_KVC2.Px()+p_KVC2.Py()*p_KVC2.Py()+p_KVC2.Pz()*p_KVC2.Pz());
      double KVC_beta = GetBeta(KVC_pdg,KVC_mom);
      
      if(KVC_beta > KVC_thre){
	
	kvc2 = true;
      }
    }

    if(mul2 || forward2){
      if(!kvc2){
	kp_count++;
      }
    }
    
  }
  double Kbeam_percent = (double)Kbeam_count / (double)tree2->GetEntries();
  cout<<Kbeam_percent<<endl;

  double kp_percent = (double)kp_count / (double)tree2->GetEntries();
  cout<<kp_percent<<endl;

  double freq = (intensity * Kbeam_percent / spill_length);
  cout<<"freq : "<<freq<<endl;

  double cs = 20; //mb
  double mb_to_cm2 = 1e-27;
  static const double Rbeam = 0.69;
  double Nreaction = intensity * kp_percent * Ntarget * mb_to_cm2*Rbeam;

  double kp_freq = Nreaction / spill_length;
  cout<<"kp freq : "<<kp_freq<<endl;
  
  
}
