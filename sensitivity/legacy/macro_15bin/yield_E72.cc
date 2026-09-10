static const int n_mombin = 15;

static const double mom_kaon[n_mombin] = {724., 726., 728., 730., 732., 734., 738., 742., 746., 750., 754., 758., 762., 766., 770.};
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
static const double T_scan = 0.5*24*3600/spill_length;  //s : integrated time of the E72 physics run
static const double T_fix = 5.5*24*3600/spill_length;
static const double eff_acc = 1.0; //efficency of accerelator
static const double lK = 0.35; //m : distance btw BAC and target center
static const double pK = 735; //MeV
static const double K_tau =  1.238*pow(10,-8);
static const double Rbeam = 0.69; //ratio of the target upstream face to the cross-sectional area of the beam profile

static const double Nbeam = FK*(T_scan + T_fix)*eff_acc*TMath::Exp(-1*lK/5.52)*Rbeam;


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

//Kp -> Lambda eta (DOI : 10.1007/BF02785525)
static const string Leta1_doi = "10.1007/BF02785525";
static const int Leta1_num = 27;
static const double Leta1_pK[Leta1_num] = {777, 806, 868, 926, 992, 1062, 1110, 1148, 1179, 1219, 1263, 1274, 1316, 1320, 1368, 1381, 1415, 1434, 1462, 1514, 1546, 1606, 1653, 1705, 1741, 1799, 1842};
static const double Leta1_pK_low[Leta1_num] = {757, 786, 838, 896, 962, 1037, 1085, 1123, 1154, 1194, 1246, 1249, 1299, 1295, 1350, 1356, 1397, 1414, 1442, 1494, 1526, 1586, 1633, 1685, 1718, 1775, 1818};
static const double Leta1_pK_high[Leta1_num] = {797, 826, 898, 956, 1022, 1087, 1135, 1173, 1204, 1244, 1280, 1299, 1333, 1345, 1386, 1406, 1433, 1454, 1482, 1534, 1566, 1626, 1673, 1725, 1764, 1823, 1866};
static const double Leta1_cs[Leta1_num] = {0.56, 0.22, 0.075, 0.0, 0.23, 0.245, 0.495, 0.37, 0.51, 0.28, 0.31, 0.39, 0.18, 0.275, 0.285, 0.29, 0.425, 0.37, 0.26, 0.27, 0.31, 0.31, 0.29, 0.3, 0.22, 0.19, 0.2};
static const double Leta1_cs_err[Leta1_num]={0.2, 0.16, 0.12, 0.12, 0.12, 0.14, 0.19, 0.13, 0.19, 0.09, 0.16, 0.08, 0.12, 0.07, 0.09, 0.07, 0.11, 0.1, 0.07, 0.08, 0.08, 0.1, 0.07, 0.08, 0.08, 0.1, 0.09};

//Kp -> Lambda eta (DOI : 10.1103/PhysRevC.64.055205)
static const string Leta2_doi = "10.1103/PhysRevC.64.055205";
static const int Leta2_num = 17;
static const double Leta2_pK[Leta2_num] = {720, 722, 724, 726, 728, 730, 732, 734, 738, 742, 746, 750, 754, 758, 762, 766, 770};
static const double Leta2_pK_err[Leta2_num] = {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5};
static const double Leta2_cs[Leta2_num] = {0.070, 0.040, 0.210, 0.470, 0.600, 0.830, 0.900, 1.150, 1.300, 1.440, 1.410, 1.240, 1.230, 1.230, 0.960, 1.020, 0.730};
static const double Leta2_cs_err[Leta2_num] = {0.020, 0.040, 0.040, 0.100, 0.070, 0.100, 0.110, 0.100, 0.100, 0.100, 0.090, 0.080, 0.070, 0.080, 0.080, 0.130, 0.190};

static const int Leta_num = Leta1_num + Leta2_num + 2;
static const double Leta_pK[Leta_num] = {700, 718, 720, 722, 724, 726, 728, 730, 732, 734, 738, 742, 746, 750, 754, 758, 762, 766, 770, 777, 806, 868, 926, 992, 1062, 1110, 1148, 1179, 1219, 1263, 1274, 1316, 1320, 1368, 1381, 1415, 1434, 1462, 1514, 1546, 1606, 1653, 1705, 1741, 1799, 1842};
static const double Leta_pK_err[Leta_num] = {2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 20.,20.,30.,30.,30.,25.,25.,25.,25.,25.,25.,25.,25.,20.};


static const double Leta_cs[Leta_num] = {0., 0., 0.070, 0.040, 0.210, 0.470, 0.600, 0.830, 0.900, 1.150, 1.300, 1.440, 1.410, 1.240, 1.230, 1.230, 0.960, 1.020, 0.730, 0.56, 0.22, 0.075, 0.0, 0.23, 0.245, 0.495, 0.37, 0.51, 0.28, 0.31, 0.39, 0.18, 0.275, 0.285, 0.29, 0.425, 0.37, 0.26, 0.27, 0.31, 0.31, 0.29, 0.3, 0.22, 0.19, 0.2};

static const double Leta_cs_err[Leta_num] = {0., 0.,0.020, 0.040, 0.040, 0.100, 0.070, 0.100, 0.110, 0.100, 0.100, 0.100, 0.090, 0.080, 0.070, 0.080, 0.080, 0.130, 0.190, 0.2, 0.16, 0.12, 0.12, 0.12, 0.14, 0.19, 0.13, 0.19, 0.09, 0.16, 0.08, 0.12, 0.07, 0.09, 0.07, 0.11, 0.1, 0.07, 0.08, 0.08, 0.1, 0.07, 0.08, 0.08, 0.1, 0.09};

double pK_to_roots(double pK){
  return sqrt((mK*mK+pK*pK)+mp*mp+2*mp*sqrt(mK*mK+pK*pK)-pK*pK);
}


static const double p_start = 600;
static const double p_end = 800;
const int np = 100.;

/*
static const double p_start = 500;
static const double p_end = 1000;
const int np = 1000.;
*/
static const double p_step = (p_end - p_start)/np;

static const double display_start = 715.;
static const double display_end = 800.;

static const double K_mean = 906;
static const double K_sigma = 13;


void yield_E72(){
  TFile *file_beam = new TFile("data/beam.k.run69_0130.root");
  

  //Beam spill calculation
  TH1D *hist_beam[6];
  TH1D *hist_beamf = new TH1D("hist_beamf","hist_beamf",np,p_start,p_end);
  TGraph *spill_graph = new TGraph(3,pK_intensity,intensity);
  
  TF1 *f_beam[6];
  //momentum scan
  for(int i=0;i<5;i++){
    hist_beam[i]= new TH1D(Form("hist_beam%d",i),Form("hist_beam%d",i),np,p_start,p_end);
    f_beam[i] = new TF1(Form("f_beam%d",i), "gaus", beam_mom[i] - 5 * K_sigma*beam_mom[i]/K_mean, beam_mom[i] + 5 * K_sigma*beam_mom[i]/K_mean);
    f_beam[i]->SetParameters(1, beam_mom[i], K_sigma*beam_mom[i]/K_mean);
    double ctaup_m = 299792458*K_tau*beam_mom[i]/mK;
    
    hist_beam[i]->FillRandom(Form("f_beam%d",i),spill_graph->Eval(beam_mom[i])*T_scan*eff_acc*TMath::Exp(-1*lK/ctaup_m)*Rbeam);
    hist_beamf->Add(hist_beam[i]);
  }
  
  //Fix beam momentum -> 735 MeV/c
  hist_beam[5] = new TH1D("hist_beam5","hist_beam5",np,p_start,p_end);
  f_beam[5] = new TF1("f_beam5","gaus",beam_mom_fix -5 * K_sigma*beam_mom_fix/K_mean, beam_mom_fix + 5 * K_sigma*beam_mom_fix/K_mean);
  f_beam[5]->SetParameters(1, beam_mom_fix, K_sigma*beam_mom_fix/K_mean);
  double ctaup_m = 299792458*K_tau*beam_mom_fix/mK;
  
  hist_beam[5]->FillRandom("f_beam5",spill_graph->Eval(beam_mom_fix)*T_fix*eff_acc*TMath::Exp(-1*lK/ctaup_m)*Rbeam);
  hist_beamf->Add(hist_beam[5]);

  TCanvas *c10 = new TCanvas("c10","c10");
  hist_beamf->Draw();

  for(int i=0;i<6;i++){
    hist_beam[i]->SetLineColor(i+1);
    hist_beam[i]->Draw("same");
  }
  /*
  
  TFile *file_nbeam = new TFile("data/n_kaon.root","READ");
  TH1D *hist_beamall = (TH1D*)file_nbeam->Get("mom_all");
  TCanvas *c10 = new TCanvas("c10","c10");
  hist_beamall->Draw();
  double binwidth = hist_beamall->GetXaxis()->GetBinWidth(1);
  cout<<"bin width : "<<binwidth<<endl;

  double xmin = hist_beamall->GetXaxis()->GetXmin();
  double xmax = hist_beamall->GetXaxis()->GetXmax();
  cout<<"xmin : "<<xmin<<", xmax : "<<xmax<<endl;
  
  */

  
  

  //Cross section data
  double Leta1_pK_err[Leta1_num];
  for(int i=0;i<Leta1_num;i++){
    Leta1_pK_err[i] = (Leta1_pK_high[i] - Leta1_pK_low[i])/2.;
  }

  
  TGraphErrors *Leta1 = new TGraphErrors(Leta1_num, Leta1_pK, Leta1_cs, Leta1_pK_err, Leta1_cs_err);
  TGraphErrors *Leta2 = new TGraphErrors(Leta2_num, Leta2_pK, Leta2_cs, Leta2_pK_err, Leta2_cs_err);
  TGraphErrors *Leta_total = new TGraphErrors(Leta_num,Leta_pK,Leta_cs,Leta_pK_err,Leta_cs_err);
  TMultiGraph *mg_cs = new TMultiGraph();
  mg_cs->Add(Leta1);
  mg_cs->Add(Leta2);

  TGraph *cs_val = new TGraph();

  for(int i=0;i<np;i++){
    double p = p_start + p_step/2. + i*p_step;
    cs_val->SetPoint(cs_val->GetN(),p,Leta_total->Eval(p));
    
  }
			   


  
  TCanvas *c1 = new TCanvas("c1","c1");
  double cs_max = 1.6;
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  c1->SetTopMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetGrid();
  Leta1->SetMarkerStyle(20);
  Leta2->SetMarkerStyle(22);
  mg_cs->SetTitle(";#it{p_{K}} [MeV/#it{c}];#sigma [mb]");
  mg_cs->GetXaxis()->SetLimits(display_start,display_end);
  mg_cs->GetYaxis()->SetRangeUser(0.,cs_max);
  
  mg_cs->Draw("AP");

  TLegend *le = new TLegend(0.8,0.5,0.48,0.6);
  le->AddEntry(Leta2,"A. Starostin #it{et al.}");
  le->AddEntry(Leta1,"R. Rader #it{et al.}");
  
  le->Draw("same");
  //show sqrt s
  TF1 *sqrt_s = new TF1("sqrt_s","x",pK_to_roots(display_start),pK_to_roots(display_end));


  TGaxis *Ecm_axis_cs = new TGaxis(display_start,cs_max,display_end,cs_max,"sqrt_s",510,"-");
  
  Ecm_axis_cs->SetTitle("#sqrt{s} [MeV]"); 
  Ecm_axis_cs->SetTitleOffset(1.4);
  Ecm_axis_cs->SetLabelFont(42);
  Ecm_axis_cs->SetTitleFont(42);
  Ecm_axis_cs->SetTitleColor(kRed);
  Ecm_axis_cs->SetLineColor(kRed);
  Ecm_axis_cs->SetLabelColor(kRed);

  Ecm_axis_cs->SetLabelOffset(0.03);
  //Ecm_axis_cs->Draw("same");
  //Yield Calculation
  double mom_cal[np];
  double mom_err[np];
  double yield_cal[np];
  double yield_err[np];
  double lumi[np];
  double lumi_err[np];

  for(int i=0;i<np;i++){
    double p = p_start + p_step/2. + i*p_step;
    mom_cal[i] = p;
    mom_err[i] = p_step/2.;
    double cs = Leta_total->Eval(p);
    double Nbeam_point = hist_beamf->GetBinContent(i+1);
    //double Nbeam_point = hist_beamall->GetBinContent(i+1);
    yield_cal[i] = Nbeam_point*Ntarget*cs*mb_to_cm2*E_eff*trigger_eff*fraction_L;
    yield_err[i] = TMath::Sqrt(yield_cal[i]);
    lumi[i] = Nbeam_point*Ntarget*1e-27*E_eff;
    cout<<lumi[i]<<endl;
    lumi_err[i] = 0.;
    

  }

  TGraphErrors *yield_total = new TGraphErrors(np,mom_cal,yield_cal,mom_err,yield_err);
  TGraphErrors *lumi_total = new TGraphErrors(np,mom_cal,lumi,mom_err,lumi_err);

  TCanvas *c5 = new TCanvas("c5","c5");
  double yield_max = 35000;
  c5->SetLeftMargin(0.15);
  c5->SetRightMargin(0.15);
  c5->SetTopMargin(0.15);
  c5->SetBottomMargin(0.15);
  c5->SetGrid();
  yield_total->SetMarkerStyle(21);
  yield_total->SetMarkerSize(0.7);
  yield_total->GetXaxis()->SetLimits(display_start,display_end);
  yield_total->GetYaxis()->SetRangeUser(0.,yield_max);
  yield_total->SetTitle(";#it{p_{K}} [MeV/#it{c}];Yield");
  yield_total->Draw("AP");
  //show sqrt s
  TGaxis *Ecm_axis_y = new TGaxis(display_start,yield_max,display_end,yield_max,"sqrt_s",510,"-");
  
  Ecm_axis_y->SetTitle("#sqrt{s} [MeV]"); 
  Ecm_axis_y->SetTitleOffset(1.4);
  Ecm_axis_y->SetLabelFont(42);
  Ecm_axis_y->SetTitleFont(42);
  Ecm_axis_y->SetTitleColor(kRed);
  Ecm_axis_y->SetLineColor(kRed);
  Ecm_axis_y->SetLabelColor(kRed);
  Ecm_axis_y->SetLabelOffset(0.03);
  //Ecm_axis_y->Draw("same");
  //Ecm_axis_y->Draw("same");

  cout<<"N_beam : "<<Nbeam<<endl;
  cout<<"N_target : "<<Ntarget<<endl;
  cout<<"Total yield : "<<yield_total->Integral()<<endl;


  ofstream outfile("param/entry_mom.txt");
  outfile << "#pK pK_min pK_max Entry" << "\n";

  ofstream outfile1("param/luminosity_mom.txt");
  outfile1 << "#pK pK_min pK_max Luminosity" << "\n";
  
  for(int i=0;i<n_mombin;i++){
    double mom_min = mom_kaon[i]-mom_kaon_err;
    double mom_max = mom_kaon[i]+mom_kaon_err;

    double entry_min = yield_total->Eval(mom_min);
    double entry_max = yield_total->Eval(mom_max);

    double lumi_min = lumi_total->Eval(mom_min);
    double lumi_max = lumi_total->Eval(mom_max);

    double entry_mom = 0.5 * (entry_min + entry_max) / (mom_max - mom_min);
    double lumi_mom = 0.5 * (lumi_min + lumi_max) / (mom_max - mom_min);
    outfile << mom_kaon[i] <<" "<<mom_min<<" "<<mom_max<<" "<<entry_mom<<"\n";
    outfile1 << mom_kaon[i] <<" "<<mom_min<<" "<<mom_max<<" "<<lumi_mom<<"\n";
    
  }

  
}
