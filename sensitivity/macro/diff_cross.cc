#include "RootData.h"
static const double xmin = -1.;
static const double xmax = 1.;
static const int n_cosbin = 15;
static const int n_mombin = 15;
static const double step = (xmax - xmin)/(double)n_cosbin;

static const double mom_kaon[n_mombin] = {724., 726., 728., 730., 732., 734., 738., 742., 746., 750., 754., 758., 762., 766., 770.};
static const double mom_err = 1.;


void diff_cross(){
  gStyle->SetOptStat(0);
 
  TFile *file = new TFile("data/mc/7201_7202_CSOn.root");
  //TFile *file = new TFile("data/mc/7201_7202_Acceptance_study.root");
  TTree *tree = (TTree*)file->Get("g4hyptpc_light");

  double mom_kaon_lab;
  double cos_theta;
  int trig_flag;
  int decay_particle_code;

  tree->SetBranchAddress("mom_kaon_lab",&mom_kaon_lab);
  tree->SetBranchAddress("cos_theta",&cos_theta);
  tree->SetBranchAddress("trig_flag",&trig_flag);
  tree->SetBranchAddress("decay_particle_code",&decay_particle_code);



  double entry[n_mombin];
  RootData::readMomentum("results/yield.root","entry_mom",n_mombin,mom_kaon,entry);

  int pass[n_mombin][n_cosbin]={0};
  int pass_mom[n_mombin]={0};
  for(int n=0;n<tree->GetEntries();n++){
    tree->GetEntry(n);
    if(decay_particle_code==2112)continue;
    if(mom_kaon_lab < mom_kaon[0] - mom_err || mom_kaon_lab > mom_kaon[n_mombin-1] + mom_err)continue;
    double cos_L = -1*cos_theta;
    int cosbin = static_cast<int>(trunc((cos_theta+1.)/step));
    int mombin = -1;
    for(int i=0;i<n_mombin;i++){
      if(mom_kaon_lab >= mom_kaon[i] - mom_err && mom_kaon_lab < mom_kaon[i] + mom_err){
	mombin = i;
	break;
      }
    }

    if(mombin < 0 || cosbin < 0 || cosbin > n_cosbin) continue;

    if(trig_flag!=0){
      pass_mom[mombin]++;
      if(cosbin==n_cosbin)pass[mombin][n_cosbin-1]++;
      else{
	pass[mombin][cosbin]++;
      }
    }
    
  }


  
  double cos_cm[n_cosbin];
  double cos_cm_err[n_cosbin];
  double acc[n_mombin][n_cosbin];
  double acc_err[n_mombin][n_cosbin];


  TGraphErrors *Leta_accept[n_mombin];
  for(int i=0;i<n_mombin;i++){
    double scale = (double)entry[i] / (double)pass_mom[i];

    for(int j=0;j<n_cosbin;j++){
      cos_cm[j] = xmin+step/2. + step*j;
      cos_cm_err[j] = step/2.;

      
      acc[i][j] = (double)pass[i][j] * scale;
      acc_err[i][j] = 0.;
      
    
    }
    Leta_accept[i] = new TGraphErrors(n_cosbin,cos_cm,acc[i],cos_cm_err,acc_err[i]);
  }

  
  
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  c1->SetTopMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetGrid();
  c1->Divide(5,3);
  for(int i=0;i<n_mombin;i++){
    c1->cd(i+1);
    Leta_accept[i]->SetTitle(";cos#theta_{#Lambda}^{CM};Total Entry");
    Leta_accept[i]->GetXaxis()->SetLimits(xmin,xmax);
    //Leta_accept[i]->GetYaxis()->SetRangeUser(0.6,0.9);
    Leta_accept[i]->SetMarkerStyle(20);
    Leta_accept[i]->SetMarkerSize(0.5);
    Leta_accept[i]->Draw("AP");
  }


  
  gSystem->mkdir("results", kTRUE);
  TFile output("results/diff_cross.root", "RECREATE");
  for (int i=0;i<n_mombin;++i) {
    RootData::momentum(output,mom_kaon[i]);
    RootData::angular(Leta_accept[i],"yield");
    TH1D passed("mc_pass","Triggered MC;cos theta;Events",n_cosbin,-1,1);
    passed.SetDirectory(nullptr);
    for (int j=0;j<n_cosbin;++j) {
      passed.SetBinContent(j+1,pass[i][j]);
      passed.SetBinError(j+1,std::sqrt(pass[i][j]));
    }
    passed.Write();
  }
  output.cd();
  output.Close();
}
