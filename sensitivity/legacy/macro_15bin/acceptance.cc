static const double xmin = -1.;
static const double xmax = 1.;
static const int n_cosbin = 15;
static const int n_mombin = 15;
static const double step = (xmax - xmin)/(double)n_cosbin;

static const double mom_kaon[n_mombin] = {724., 726., 728., 730., 732., 734., 738., 742., 746., 750., 754., 758., 762., 766., 770.};
static const double mom_err = 1.;


void acceptance(){
  gStyle->SetOptStat(0);
  
  TFile *file = new TFile("../data/7201_7202_Acceptance_study.root");
  //TFile *file = new TFile("../data/mom735_flat_eta_lambda.root");
  //TFile *file = new TFile("../data/7201_7202_CSOn.root");
  TTree *tree = (TTree*)file->Get("g4hyptpc_light");

  double mom_kaon_lab;
  double cos_theta;
  int trig_flag;
  int decay_particle_code;

  tree->SetBranchAddress("mom_kaon_lab",&mom_kaon_lab);
  tree->SetBranchAddress("cos_theta",&cos_theta);
  tree->SetBranchAddress("trig_flag",&trig_flag);
  tree->SetBranchAddress("decay_particle_code",&decay_particle_code);

  
  int pass[n_mombin][n_cosbin]={0};
  int total[n_mombin][n_cosbin]={0};

  int pass_all[n_cosbin]={0};
  int total_all[n_cosbin]={0};
  for(int n=0;n<tree->GetEntries();n++){
    tree->GetEntry(n);
    if(decay_particle_code==2112)continue;

    int cosbin = static_cast<int>(trunc((cos_theta+1.)/step));
    
    if(cosbin==n_cosbin)total_all[n_cosbin-1]++;
    else{total_all[cosbin]++;}

    if(trig_flag!=0){
      if(cosbin==n_cosbin)pass_all[n_cosbin-1]++;
      else{pass_all[cosbin]++;}
    }
    
    if(mom_kaon_lab < mom_kaon[0] - mom_err || mom_kaon_lab > mom_kaon[n_mombin-1] + mom_err)continue;

    int mombin;
    for(int i=0;i<n_mombin;i++){
      if(mom_kaon_lab >= mom_kaon[i] - mom_err && mom_kaon_lab < mom_kaon[i] + mom_err){
	mombin = i;
	break;
      }
    }
    
    if(cosbin==n_cosbin)total[mombin][n_cosbin-1]++;
    else{total[mombin][cosbin]++;}

    if(trig_flag!=0){
      if(cosbin==n_cosbin)pass[mombin][n_cosbin-1]++;
      else{pass[mombin][cosbin]++;}
    }
    
  }
  double cos_cm[n_cosbin];
  double cos_cm_err[n_cosbin];
  double acc[n_mombin][n_cosbin];
  double acc_err[n_mombin][n_cosbin];
  double acc_all[n_cosbin];
  double acc_all_err[n_cosbin];


  TGraphErrors *Leta_accept[n_mombin];
  TGraphErrors *Leta_accept_all;


  ofstream outfile("param/acceptance_mom_cos.txt");
  outfile <<"#pK pK_min pK_max cos_cm cos_cm_err Acceptance"<< "\n";
  for(int i=0;i<n_mombin;i++){
    for(int j=0;j<n_cosbin;j++){
      cos_cm[j] = xmin+step/2. + step*j;
      cos_cm_err[j] = step/2.;

      acc[i][j] = (double)pass[i][j]/(double)total[i][j];
      acc_err[i][j] = TMath::Sqrt(acc[i][j]*0.01*(1.-acc[i][j]*0.01)/(double)total[i][j]);
      outfile << mom_kaon[i] <<" "<<mom_kaon[i]-mom_err<<" "<<mom_kaon[i]+mom_err<<" "<<cos_cm[j]<<" "<<cos_cm_err[j]<<" "<<acc[i][j]<<"\n";
      if(i==0){
	acc_all[j] = (double)pass_all[j]/(double)total_all[j];
	acc_all_err[j] = TMath::Sqrt(acc_all[j]*0.01*(1.-acc_all[j]*0.01)/(double)total_all[j]);
	
      }
	
    
    }
    Leta_accept[i] = new TGraphErrors(n_cosbin,cos_cm,acc[i],cos_cm_err,acc_err[i]);
    if(i==0) Leta_accept_all = new TGraphErrors(n_cosbin,cos_cm,acc_all,cos_cm_err,acc_all_err);
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
    Leta_accept[i]->SetTitle(";cos#theta_{#eta}^{CM};Acceptance");
    Leta_accept[i]->GetXaxis()->SetLimits(xmin,xmax);
    Leta_accept[i]->GetYaxis()->SetRangeUser(0.6,1.);
    Leta_accept[i]->SetMarkerStyle(20);
    Leta_accept[i]->SetMarkerSize(0.5);
    Leta_accept[i]->Draw("AP");
  }


  TCanvas *c2= new TCanvas("c2","c2");
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.15);
  c2->SetTopMargin(0.15);
  c2->SetBottomMargin(0.15);
  gPad->SetGrid();
  Leta_accept_all->SetTitle(";cos#theta_{#eta}^{CM};Acceptance");
  Leta_accept_all->GetXaxis()->SetLimits(xmin,xmax);
  Leta_accept_all->GetYaxis()->SetRangeUser(0.6,1.);
  Leta_accept_all->SetMarkerStyle(20);
  Leta_accept_all->SetMarkerSize(1);
  Leta_accept_all->Draw("AP");
  
  
  
}
