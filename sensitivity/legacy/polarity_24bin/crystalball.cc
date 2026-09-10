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
static const double K_tau =  1.238*pow(10,-8);
static const double Rbeam = 0.69; //ratio of the target upstream face to the cross-sectional area of the beam profile

static const double Nbeam = FK*(T_scan + T_fix)*eff_acc*TMath::Exp(-1*lK/5.52)*Rbeam;


static const double rho = 0.07085; //g/cm3

static const double NA = 6.022 *1e23; //mol-1
static const double W = 1.008; //atomic mass of LH2 target
static const double Ltarget = 7.139; //cm targe thickness-> approximately assume

static const double Ntarget = rho*NA/W*Ltarget;

static const double E_eff = 0.8; //DAQ effciency * offline analysis efficiency
static const double trigger_eff = 0.9; //HToF - MP1?
static const double fraction_L = 0.64; //Lambda -> ppi
static const double mb_to_cm2 = 1e-27;


static const double xmin = -1.;
static const double xmax = 1.;
static const int n_cosbin = 15;
//static const int n_mombin = 15;
static const int n_mombin = 24;
static const double step = (xmax - xmin)/(double)n_cosbin;

static const int numpK = 15;
static const int numcosTheta = 9;
static const double pK[numpK] = { 724, 726, 728, 730, 732, 734, 738, 742, 746, 750, 754, 758, 762, 766, 770 };
static const double pK_err[numpK] = { 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5 };
static const double cosTheta[numcosTheta] = { -0.95, -0.80, -0.60, -0.35, 0.00, 0.35, 0.60, 0.80, 0.95 };
static const double cosTheta_err[numcosTheta] = { 0.05, 0.10, 0.10, 0.15, 0.20, 0.15, 0.10, 0.10, 0.05 };
static const double diff_cs[numpK][numcosTheta] = {{ 0.01, 0.03, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03 },
						   { 0.07, 0.02, 0.04, 0.05, 0.02, 0.02, 0.06, 0.05, 0.03 },
						   { 0.07, 0.05, 0.04, 0.05, 0.04, 0.02, 0.08, 0.08, 0.02 },
						   { 0.13, 0.08, 0.04, 0.04, 0.05, 0.06, 0.07, 0.10, 0.14 },
						   { 0.12, 0.08, 0.11, 0.04, 0.06, 0.06, 0.06, 0.09, 0.08 },
						   { 0.17, 0.11, 0.10, 0.04, 0.06, 0.07, 0.09, 0.16, 0.23 },
						   { 0.13, 0.14, 0.08, 0.07, 0.07, 0.12, 0.13, 0.12, 0.17 },
						   { 0.16, 0.16, 0.11, 0.11, 0.08, 0.10, 0.15, 0.11, 0.15 },
						   { 0.17, 0.13, 0.14, 0.09, 0.11, 0.09, 0.11, 0.13, 0.09 },
						   { 0.14, 0.09, 0.12, 0.09, 0.07, 0.08, 0.10, 0.11, 0.14 },
						   { 0.15, 0.12, 0.09, 0.10, 0.10, 0.09, 0.08, 0.11, 0.10 },
						   { 0.13, 0.13, 0.10, 0.11, 0.06, 0.09, 0.11, 0.10, 0.11 },
						   { 0.11, 0.11, 0.10, 0.09, 0.07, 0.05, 0.08, 0.06, 0.02 },
						   { 0.11, 0.12, 0.14, 0.09, 0.08, 0.05, 0.04, 0.08, 0.04 },
						   { 0.06, 0.06, 0.11, 0.03, 0.08, 0.03, 0.07, 0.11, 0.06 }};

static const double diff_cs_err[numpK][numcosTheta] = {{ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02 },
						       { 0.05, 0.02, 0.03, 0.02, 0.01, 0.01, 0.03, 0.03, 0.03 },
						       { 0.03, 0.02, 0.02, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02 },
						       { 0.05, 0.03, 0.02, 0.02, 0.02, 0.02, 0.03, 0.03, 0.05 },
						       { 0.04, 0.03, 0.04, 0.02, 0.02, 0.02, 0.03, 0.03, 0.04 },
						       { 0.05, 0.03, 0.03, 0.01, 0.02, 0.02, 0.03, 0.03, 0.06 },
						       { 0.04, 0.03, 0.02, 0.02, 0.02, 0.02, 0.03, 0.03, 0.04 },
						       { 0.04, 0.03, 0.02, 0.02, 0.01, 0.02, 0.03, 0.02, 0.04 },
						       { 0.04, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.03 },
						       { 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.02, 0.02, 0.03 },
						       { 0.03, 0.02, 0.02, 0.02, 0.01, 0.01, 0.02, 0.02, 0.02 },
						       { 0.03, 0.02, 0.02, 0.02, 0.01, 0.02, 0.02, 0.02, 0.03 },
						       { 0.03, 0.02, 0.02, 0.02, 0.02, 0.01, 0.02, 0.02, 0.01 },
						       { 0.05, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.03, 0.03 },
						       { 0.06, 0.04, 0.06, 0.03, 0.04, 0.04, 0.05, 0.06, 0.06 }};

void convert_roots(double M1, double P1, double M2, double P2, double M1_err, double P1_err, double M2_err, double P2_err, double &roots, double &roots_err){
  double E1 = TMath::Sqrt(M1*M1 + P1*P1);
  double E2 = TMath::Sqrt(M2*M2 + P2*P2);
  double totalE = E1 + E2;
  double totalpz = P1 + P2;

  double error_E1 = TMath::Sqrt((M1 * M1 * M1_err * M1_err) / (4 * E1 * E1) +
                                   (P1 * P1 * P1_err * P1_err) / (4 * E1 * E1));
  double error_E2 = TMath::Sqrt((M2 * M2 * M2_err * M2_err) / (4 * E2 * E2) +
                                   (P2 * P2 * P2_err * P2_err) / (4 * E2 * E2));

  double error_totalE = TMath::Sqrt(error_E1 * error_E1 +  error_E2 * error_E2);
  double error_totalpz = TMath::Sqrt(P1_err * P1_err + P2_err * P2_err);

  roots =  TMath::Sqrt(totalE*totalE - totalpz*totalpz);

  /*
  roots_err =TMath::Sqrt((totalE * totalE * error_totalE * error_totalE) +
                                           (totalpz * totalpz * error_totalpz * error_totalpz)) /( totalE);
  */

  roots_err=0;


}

double roots_mom(double M1, double M2, double roots){
  return TMath::Sqrt(pow((roots*roots-M1*M1-M2*M2),2)/(4*M2*M2)-M1*M1);
  
}

double Legendre2(double *x, double *par){
  return par[0] + par[1]*x[0] + par[2]*(3/2*x[0]*x[0]-1/2);
}

double Legendre0(double *x, double *par){
  return par[0];
}


double Legendre2(double x, double *par){
  return par[0] + par[1]*x + par[2]*(3/2*x*x-1/2);
}

double fit_amp(double *x1, double *par1){
  return par1[1]/(2*TMath::Pi()*pow(x1[0]-par1[0],2)+pow(par1[1]/2,2));
}


double GlobalFitFunction(double* x, double* par) {
    int idx = static_cast<int>(x[1]);  // Index to differentiate between A0, A1, A2 (0 for A0, 1 for A1, 2 for A2)
    double E = x[0];      // Energy
    double M = par[0];    // Resonance mass (공통)
    double Gamma = par[1]; // Resonance width (공통)

    // Scale factors for A0, A1, A2
    double A0_scale = par[2];
    double A1_scale = par[3];
    double A2_scale = par[4];
    
    // Breit-Wigner formula
    double commonBW = (Gamma*Gamma) / ((E - M)*(E - M) + (Gamma*Gamma)/4.0);
    
    // Return different A_i based on idx
    if (idx == 0) return A0_scale * commonBW;  // A0
    if (idx == 1) return A1_scale * commonBW;  // A1
    return A2_scale * commonBW;  // A2
}

void crystalball(){

  double roots[numpK];
  double roots_err[numpK];
  for(int i=0;i<numpK;i++){
    convert_roots(mp, 0, mK, pK[i], 0, 0, 0, pK_err[i], roots[i], roots_err[i]);

  }

  

  TGraphErrors *diff_graph[numpK];
  TF1 *poly0[numpK];
  TF1 *poly[numpK];
  double par_legen0[numpK][3];
  double par_legen_err0[numpK][3];
  double par_legen[numpK][3];
  double par_legen_err[numpK][3];
  for(int i=0;i<numpK;i++){
    diff_graph[i] = new TGraphErrors(numcosTheta, cosTheta, diff_cs[i], cosTheta_err, diff_cs_err[i]);
    diff_graph[i]->GetXaxis()->SetRangeUser(-1,1);
    diff_graph[i]->GetYaxis()->SetRangeUser(0,0.3);
    diff_graph[i]->SetMarkerStyle(21);
    diff_graph[i]->SetMarkerSize(0.5);
    diff_graph[i]->GetXaxis()->SetTitleSize(0.08);
    diff_graph[i]->GetXaxis()->SetLabelSize(0.08);
    diff_graph[i]->GetXaxis()->SetTitleOffset(0.8);
    diff_graph[i]->GetYaxis()->SetTitleSize(0.06);
    diff_graph[i]->GetYaxis()->SetLabelSize(0.08);
    diff_graph[i]->GetYaxis()->SetTitleOffset(0.8);
    diff_graph[i]->GetYaxis()->SetNdivisions(405);

    

    poly0[i] = new TF1(Form("poly0%d",i),Legendre0,-1,1,3);
    poly0[i]->SetLineColor(kBlue);
    poly0[i]->SetLineWidth(3);
    
    poly[i] = new TF1(Form("poly%d",i),Legendre2,-1,1,3);
    poly[i]->SetLineColor(kRed);
    poly[i]->SetLineWidth(3);


    diff_graph[i]->Fit(poly[i],"Q0","",-1,1);
    diff_graph[i]->Fit(poly0[i],"Q0","",-1,1);
    
    poly0[i]->GetParameters(par_legen0[i]);
    poly[i]->GetParameters(par_legen[i]);

    for(int j=0;j<3;j++){
      par_legen_err0[i][j] = poly0[i]->GetParError(j);
      par_legen_err[i][j] = poly[i]->GetParError(j);
    }
  }

  double par_legen_swap0[3][numpK];
  double par_legen_swap_err0[3][numpK];
  
  double par_legen_swap[3][numpK];
  double par_legen_swap_err[3][numpK];
  for (int i = 0; i < numpK; ++i) {
    for (int j = 0; j < 3; ++j) {
      par_legen_swap0[j][i] = par_legen0[i][j];
      par_legen_swap_err0[j][i] = par_legen_err0[i][j];
      
      par_legen_swap[j][i] = par_legen[i][j];
      par_legen_swap_err[j][i] = par_legen_err[i][j];
    }
  }

  TGraphErrors *Amp[3];
  TGraphErrors *Amp0[3];
  TF1 *f_amp[3];
  
  for(int i=0;i<3;i++){
    Amp0[i] = new TGraphErrors(numpK, roots, par_legen_swap0[i], roots_err, par_legen_swap_err0[i]);
    Amp[i] = new TGraphErrors(numpK, roots, par_legen_swap[i], roots_err, par_legen_swap_err[i]);
    f_amp[i] = new TF1(Form("f_amp%d",i),fit_amp,1660,1690,2);
    
  }

  

  TCanvas *c1 = new TCanvas("c1","c1",269,53,1369,540);
  //gStyle->SetTitleFontSize(0.2);
  c1->Divide(5,3);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  c1->SetBottomMargin(0.15);
  for(int i=0;i<numpK;i++){
    c1->cd(i+1);
    //diff_graph[i]->SetTitle(Form("W = %.1f;cos #theta;d#sigma/d#Omega [mb/sr]",roots[i]));
    diff_graph[i]->SetTitle(Form("W = %.1f;;",roots[i]));
    


    diff_graph[i]->Draw("AP");
    poly[i]->Draw("same");
    //poly0[i]->Draw("same");
  }

  
  TCanvas *c2[3];
  for(int i=0;i<3;i++){
    c2[i]= new TCanvas(Form("c2%d",i),Form("c2%d",i),800,650);

    c2[i]->SetTopMargin(0.15);
    c2[i]->SetBottomMargin(0.15);
    c2[i]->SetRightMargin(0.15);
    c2[i]->SetLeftMargin(0.15);
  }
  /*
  auto mg = new TMultiGraph();
  Amp0[0]->SetMarkerStyle(20);
  Amp0[0]->SetMarkerSize(1);
  Amp[0]->SetMarkerStyle(21);
  Amp[0]->SetMarkerSize(1);
  mg->Add(Amp0[0]);
  mg->Add(Amp[0]);
  */


  //TF1* globalFit = new TF1("globalFit", GlobalFitFunction, 1660., 1690., 5);
  TF1* globalFit = new TF1("globalFit", GlobalFitFunction, 1660., 1680., 5); 
    
    // 파라미터 초기값 설정
  globalFit->SetParameters(1670.8, 7.9, 0.105, 0.85, 0.026);
  globalFit->SetParLimits(0,1668,1675);
  globalFit->SetParLimits(1,5,10);
  globalFit->SetParLimits(2,0.,0.15);
  globalFit->SetParLimits(3,0.5,1.);
  globalFit->SetParLimits(4,0.,0.1);
  

  /*
    TVirtualFitter::Fitter(0, 5);  // 5개의 파라미터 사용
    globalFit->SetParNames("M", "Gamma", "A0_scale", "A1_scale", "A2_scale");

    // Fit A0, A1, A2 동시에
    globalFit->SetNpx(numpK*10); // 정확도 증가를 위한 세부 조정
    Amp[0]->Fit("globalFit", "Q"); // Fit A0 데이터
    Amp[1]->Fit("globalFit", "Q+"); // Fit A1 데이터
    Amp[2]->Fit("globalFit", "Q+"); // Fit A2 데이터

    cout<<globalFit->GetParameter(0)<<endl;
    cout<<globalFit->GetParameter(1)<<endl;
    cout<<globalFit->GetParameter(2)<<endl;
    cout<<globalFit->GetParameter(3)<<endl;
    cout<<globalFit->GetParameter(4)<<endl;
  
  for(int i=0;i<3;i++){
    c2[i]->cd();
    Amp[i]->SetTitle(Form(";#sqrt{s} [MeV];A_{%d}",i));
    Amp[i]->SetLineWidth(2);
    
    Amp[i]->GetYaxis()->SetRangeUser(-0.1,0.25);
    Amp[i]->SetMarkerStyle(21);
    Amp[i]->SetMarkerSize(1.2);

    Amp[i]->GetXaxis()->SetTitleSize(0.06);    // Set title size (default is 0.05)
    Amp[i]->GetXaxis()->SetLabelSize(0.05);    // Set label size (default is 0.04)
    Amp[i]->GetXaxis()->SetTitleOffset(1.0);   // Set title offset (default is 1)
    
    // Customize Y axis
    Amp[i]->GetYaxis()->SetTitleSize(0.06);    // Set title size
    Amp[i]->GetYaxis()->SetLabelSize(0.05);    // Set label size

    Amp[i]->GetYaxis()->SetTitleOffset(1.2);
    Amp[i]->Draw("AP");
    //    Amp[i]->Fit(f_amp[i],"","",1660,1680);
  }
  */

  

  TGraph2D * CB2D = new TGraph2D();
  for(int i=0;i<numpK;i++){
    for(int j=0;j<numcosTheta;j++){
      CB2D->AddPoint(pK[i],cosTheta[j],diff_cs[i][j]);
    }
  }
  TH2D *hist_CB2D = CB2D -> GetHistogram();

  double minpK = pK[0];
  double maxpK = pK[14];
  double mincos = -1;
  double maxcos = 1;
  
  TGraph2D *CB2D_test = new TGraph2D();
  for(int n=0;n<100000;n++){
    double randpK = minpK + (maxpK - minpK) *gRandom->Rndm();
    double randcos = mincos + (maxcos - mincos) *gRandom->Rndm();
    CB2D_test->AddPoint(randpK,randcos,CB2D->Interpolate(randpK,randcos));
  }


  TCanvas *c3 = new TCanvas("c3","c3",800,650);
  c3->Divide(3);
  c3->cd(1);
  //CB2D->Draw("surf1");
  CB2D->Draw("colz");
  c3->cd(2);
  hist_CB2D->Draw("colz");
  c3->cd(3);
  CB2D_test->Draw("colz");



  //After doing Legendre fitting
  TGraph2D *CB_fit = new TGraph2D();
  int numcostest = 10;
  for(int i=0;i<numpK;i++){
    for(int j=0;j<numcostest+1;j++){
      double pointcos = mincos + (maxcos - mincos)/numcostest*j;
      CB_fit->AddPoint(pK[i], pointcos, Legendre2(pointcos,par_legen[i]));
    }
  }


  TGraph2D *CB2D_test2 = new TGraph2D();
  for(int n=0;n<100000;n++){
    double randpK = minpK + (maxpK - minpK) *gRandom->Rndm();
    double randcos = mincos + (maxcos - mincos) *gRandom->Rndm();
    CB2D_test2->AddPoint(randpK,randcos,CB_fit->Interpolate(randpK,randcos));
  }

  TCanvas *c4 = new TCanvas("c4","c4",800,650);
  c4->Divide(3);
  c4->cd(1);
  CB2D->Draw("colz");
  c4->cd(2);
  CB_fit->Draw("colz");
  c4->cd(3);
  CB2D_test2->Draw("colz");

    
  cout<<roots_mom(mK,mp,(meta+mLambda))<<endl;
  //cout<<roots_mom(mK,mp,1670)<<endl;


    //Calculate differential Cross Section Error
  ifstream infile("param/luminosity_mom.txt");
  double luminosity[n_mombin];
  string line;
  int count = 0;
  while(getline(infile,line)){
    if(line.empty() || line[0] == '#')continue;

    istringstream iss(line);
    double mom, mom_min, mom_max, m_lumi;
    if(!(iss >> mom >> mom_min >> mom_max >> m_lumi)){
      continue;
    }
    luminosity[count] = m_lumi;
    count++;

  }

  infile.close();

  ifstream infile1("param/entry_mom_cos.txt");
  double N_mom_cos[n_mombin][n_cosbin];
  string line1;
  int count_mom = 0;
  int count_cos = 0;
  while(getline(infile1,line1)){
    if(line1.empty() || line1[0] == '#')continue;

    istringstream iss(line1);
    double mom, mom_min, mom_max, cos, cos_err, m_N;
    if(!(iss >> mom >> mom_min >> mom_max >> cos >> cos_err >> m_N)){
      continue;
    }

    N_mom_cos[count_mom][count_cos] = m_N;
    count_cos++;
    if(count_cos % n_cosbin == 0){
      count_cos = 0;
      count_mom ++;
    }
  }

  ifstream infile2("param/acceptance_mom_cos.txt");
  double acc_mom_cos[n_mombin][n_cosbin];
  string line2;
  int count2_mom = 0;
  int count2_cos = 0;
  while(getline(infile2,line2)){
    if(line2.empty() || line2[0] == '#')continue;

    istringstream iss(line2);
    double mom, mom_min, mom_max, cos, cos_err, m_acc;
    if(!(iss >> mom >> mom_min >> mom_max >> cos >> cos_err >> m_acc)){
      continue;
    }

    acc_mom_cos[count2_mom][count2_cos] = m_acc;
    count2_cos++;
    if(count2_cos % n_cosbin == 0){
      count2_cos = 0;
      count2_mom ++;
    }
  }



  
  //Experiment Expectation
  TGraphErrors* g_diff_exp[n_cosbin];
  double cos_cm[n_cosbin];
  double cos_cm_err[n_cosbin];
  double diff_cs_exp[n_mombin][n_cosbin];
  double diff_cs_exp_err[n_mombin][n_cosbin];
  for(int i=0;i<n_mombin;i++){
    for(int j=0;j<n_cosbin;j++){
      cos_cm[j] = xmin+step/2. + step*j;
      cos_cm_err[j] = step/2.;
      diff_cs_exp[i][j] = poly[i]->Eval(cos_cm[j]);
      diff_cs_exp_err[i][j] = TMath::Sqrt(N_mom_cos[i][j]) / (luminosity[i] * step * 2.0 * TMath::Pi() * acc_mom_cos[i][j]);
    }
    g_diff_exp[i] = new TGraphErrors(n_cosbin,cos_cm,diff_cs_exp[i],cos_cm_err,diff_cs_exp_err[i]);
  }

  cout<<"1"<<endl;
  
  TCanvas* c0 = new TCanvas("c0","c0",269,53,1369,540);
  c0->Divide(5,3);
  c0->SetLeftMargin(0.15);
  c0->SetRightMargin(0.15);
  c0->SetBottomMargin(0.15);
  c0->SetGrid();
  TMultiGraph *mg_diff[n_mombin];
  for(int i=0;i<n_mombin;i++){
    mg_diff[i] = new TMultiGraph();
    c0->cd(i+1);
    gPad->SetGrid();
    g_diff_exp[i]->SetMarkerStyle(20);
    g_diff_exp[i]->SetMarkerSize(0.3);
    g_diff_exp[i]->SetMarkerColor(kRed);
    g_diff_exp[i]->SetLineColor(kRed);
    
    mg_diff[i]->GetXaxis()->SetLimits(xmin,xmax);
    mg_diff[i]->GetYaxis()->SetRangeUser(0.,0.2);
    mg_diff[i]->Add(diff_graph[i]);
    mg_diff[i]->Add(g_diff_exp[i]);
    
    /*
    mg_diff[i]->GetXaxis()->SetTitleSize(0.06);
    mg_diff[i]->GetXaxis()->SetLabelSize(0.08);
    mg_diff[i]->GetXaxis()->SetTitleOffset(0.8);
    mg_diff[i]->GetYaxis()->SetTitleSize(0.08);
    mg_diff[i]->GetYaxis()->SetLabelSize(0.08);
    mg_diff[i]->GetYaxis()->SetTitleOffset(0.8);
    mg_diff[i]->GetYaxis()->SetNdivisions(405);
    mg_diff[i]->GetXaxis()->SetNdivisions(405);
    

    */
    

    mg_diff[i]->Draw("APZ");
    TLatex* latex = new TLatex();
    latex->SetNDC();  
    latex->SetTextSize(0.08);
    latex->SetTextAlign(22); 
    latex->DrawLatex(0.5, 0.95, Form("%.1f MeV",roots[i]));
  }

  


  TCanvas *c10 = new TCanvas("c10","c10");
  c10->SetLeftMargin(0.15);
  c10->SetRightMargin(0.15);
  c10->SetBottomMargin(0.15);
  c10->SetGrid();
  g_diff_exp[6]->SetMarkerStyle(21);
  g_diff_exp[6]->SetMarkerSize(0.8);
  //mg_diff[6]->Draw("APZ");
  mg_diff[6]->GetXaxis()->SetTitle("cos#theta_{#eta}^{CM}");
  mg_diff[6]->GetYaxis()->SetTitle("d#sigma/d#Omega [mb/sr]");

  mg_diff[6]->Draw("AP");


  
  
}
