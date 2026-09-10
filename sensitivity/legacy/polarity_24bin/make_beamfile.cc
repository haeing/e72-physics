static const int n_mombin = 24;
static const double mom_min = 724.;
static const double mom_max = 770.;
static const double mom_step = 2.;

double mom_kaon[n_mombin];
static const double mom_kaon_err = 1.;

static const double mK = 493.677;          //usbar, ubars
static const double mK_err = 0.016;
static const double mp = 938.27208816;     
static const double mp_err = 0.00000029;
static const double mn = 939.5654205;
static const double meta =547.862;         //c1(uubar + ddbar) + c2(ssbar)
static const double mLambda = 1115.683;    //uds
static const double mK0 = 497.611;         //dbars
static const double mSigmap = 1189.37;     //uus
static const double mSigma0 = 1192.642;    //uds
static const double mSigmam = 1197.449;    //dds
static const double mpi = 139.57039;       //udbar, dubar
static const double mpi0 = 134.9768;       //(uubar - ddbar)/sqrt(2)x

static const double pK_intensity[3] = {700,735,750};
static const double intensity[3] = {27400,38400,44000};


static const double beam_mom[5] = {685.,705.,725.,745.,765.};
static const double beam_mom_fix = 735.;

static const double FK = 3.8*pow(10,4); //# of beam particles per spill 
static const double spill_length = 4.24; // [s]
//static const double T_scan = 0.5*24*3600/spill_length;  //s : integrated time of the E72 physics run
//static const double T_fix = 5.5*24*3600/spill_length;
static const double eff_acc = 1.0; //efficency of accerelator
static const double lK = 0.35; //m : distance btw BAC and target center
static const double pK = 735; //MeV
static const double K_tau =  1.238*pow(10,-8);
static const double Rbeam = 0.69; //ratio of the target upstream face to the cross-sectional area of the beam profile
static const double beam_power_factor = 90./82.;

//static const double Nbeam = FK*(T_scan + T_fix)*eff_acc*TMath::Exp(-1*lK/5.52)*Rbeam;


static const double rho = 0.07085; //g/cm3

static const double NA = 6.022 *1e23; //mol-1
static const double W = 1.008; //atomic mass of LH2 target
static const double Ltarget = 7.139; //cm targe thickness-> approximately assume

static const double Ntarget = rho*NA/W*Ltarget;

static const double E_eff = 0.8; //DAQ effciency * offline analysis efficiency
//static const double trigger_eff = 0.9; //HTOF,BAC,KVC
static const double trigger_eff = 0.78; //HTOF,BAC,KVC 
static const double fraction_L = 0.64; //Lambda -> ppi
static const double mb_to_cm2 = 1e-27;


static const double p_start = 620;
static const double p_end = 950;
//static const double p_step = 0.5;
static const double p_step = 2;
const int np = (p_end - p_start)/p_step;
//const int np = 1000.;
//static const double p_step = (p_end - p_start)/np;

static const double display_start = 715.;
static const double display_end = 800.;

static const double K_mean = 906;
static const double K_sigma = 13;

void make_beamfile(){
  cout<<np<<endl;
  TH1D *hist_beam[5];
  for(int i=0;i<5;i++){
    hist_beam[i] = new TH1D(Form("hist_beam_%d",35+5*i),Form("hist_beam_%d",35+5*i),np,p_start,p_end);
  }
  TH1D *hist_original = new TH1D("hist_original","hist_original",np,p_start,p_end);
  TFile *file = new TFile("../data/beam_profile_run344_layer13_center.root");
  TTree *tree= (TTree*)file->Get("tr");

  TFile *file_beam = new TFile("../data/n_kaon_scan.root","recreate");
  

  double pInx;
  double pIny;
  double pInz;

  tree->SetBranchAddress("pInx",&pInx);
  tree->SetBranchAddress("pIny",&pIny);
  tree->SetBranchAddress("pInz",&pInz);

  TGraph *spill_graph = new TGraph(3,pK_intensity,intensity);
  
  for(int n=0;n<tree->GetEntries();n++){
    tree->GetEntry(n);
    hist_original->Fill(TMath::Sqrt(pInx*pInx+pIny*pIny+pInz*pInz)*1000);
  }

  cout<<hist_original->Integral()<<endl;
  for(int i=1;i<hist_original->GetNbinsX()+1;i++){
    double mom = hist_original->GetBinCenter(i);
    double entry;
    if(mom > 800)entry = 0;
    else{ entry = hist_original->GetBinContent(i);}
    //if(entry == 0)continue;
    double total_entry = hist_original->Integral();

    double ctaup_m = 299792458*K_tau*mom/mK;
    
    for(int j=0;j<5;j++){
      double T_fix = (3.5+0.5*j)*24*3600/spill_length;
      double T_scan = 7.*24*3600/spill_length - T_fix;
      double entry_per_bin_scan = (spill_graph->Eval(beam_mom_fix)*T_scan*eff_acc*TMath::Exp(-1*lK/ctaup_m)*Rbeam*beam_power_factor)/(double)np;
      hist_beam[j]->SetBinContent(i,entry/total_entry *(spill_graph->Eval(beam_mom_fix)*T_fix*eff_acc*TMath::Exp(-1*lK/ctaup_m)*Rbeam)*beam_power_factor +entry_per_bin_scan);
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(3,2);
  for(int i=0;i<5;i++){
    c1->cd(i+1);
    hist_beam[i]->Draw();
    hist_beam[i]->Write();
    cout<<3.5 + 0.5*i<<" : "<<hist_beam[i]->Integral()<<endl;
  }
  file_beam->Close();
  
  
}
