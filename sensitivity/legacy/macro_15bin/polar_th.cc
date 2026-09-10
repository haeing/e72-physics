static const double xmin = -1.;
static const double xmax = 1.;

static const int n_cosbin = 15;
static const int n_mombin = 15;
//static const int n_mombin = 24;
static const double step = (xmax - xmin)/(double)n_cosbin;

const int n_cos = 9;

const double alpha = 0.65; //Asymmetry parameter of Lambda

static const double mom_kaon[n_mombin] = {724., 726., 728., 730., 732., 734., 738., 742., 746., 750., 754., 758., 762., 766., 770.};
//static const double mom_kaon[n_mombin] = {724., 726., 728., 730., 732., 734., 736., 738., 740., 742., 744., 746., 748., 750., 752., 754., 756., 758., 760., 762., 764., 766., 768., 770.};
static const double mom_err = 1.;

void polar_th() {
  double coscm[n_cos]     = {-0.95, -0.80, -0.60, -0.35, 0.00, 0.35, 0.60, 0.80, 0.85};
  double coscm_err[n_cos] = {0.05,  0.10,  0.10,  0.15, 0.20, 0.15, 0.10, 0.10, 0.05};

  Double_t pol[2][n_cos] = {
    { 0.16, 0.27, -0.26, -0.17, 0.05, 0.05, 0.21, 0.16, -0.02},
    { 0.06, -0.16, 0.14, -0.02, -0.04, 0.59, 0.07, -0.48, -0.13}
  };

  Double_t pol_err[2][n_cos] = {
    {0.29, 0.23, 0.27, 0.22, 0.22, 0.24, 0.23, 0.23, 0.27},
    {0.32, 0.23, 0.26, 0.23, 0.24, 0.26, 0.26, 0.26, 0.38}
  };

  std::vector<double> cos_theta[2][2];
  std::vector<double> polarity[2][2];

  // === 735 MeV/c ===
  //Dwave
  cos_theta[0][0] = {-1.0, -0.4797, -0.21805, 0.1594, 0.48271, 0.73383, 1.0};
  polarity[0][0] = {0.0, -0.22509, -0.14823, 0.20741, 0.3371, 0.27142, 0.0};
  
  cos_theta[0][1] = {-1.0,-0.55779, 0.10178, 1.0};
  polarity[0][1] = {0.0, -0.64754, -0.56467, 0.00};

  // === 765 MeV/c ===
  //Dwave
  cos_theta[1][0] ={-1.0, -0.80968, -0.48915, -0.05509, 0.28548, 0.62604, 0.74624, 0.84975, 0.95659, 1.00};
  polarity[1][0] = {-0.00, -0.03167, -0.0362, 0.02715, 0.09955, 0.14932,0.1448, 0.1267, 0.0724, 0.0};
    
  cos_theta[1][1] ={-1.0, -0.43239, -0.00167, 0.28214, 0.69616, 1.0};
  polarity[1][1] = {0.0, -0.38017, -0.53157, -0.47747, -0.23689, 0.0};



  const char* mom_label[2] = {"735 MeV/c", "765 MeV/c"};
  int colors[2] = {kRed+1, kBlue+1};

  //Experiment expectation
  
  ifstream infile1("param/entry_mom_cos.txt");
  double N_mom_cos[n_mombin][n_cosbin];
  double N_mom[n_mombin] = {0.};
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
    N_mom[count_mom] += m_N;
    count_cos++;
    if(count_cos % n_cosbin == 0){
      count_cos = 0;
      count_mom ++;
    }
  }
  infile1.close();

  double N_mom_cos_compare[n_cosbin] = {0.};
  for(int i=0;i<10;i++){
    for(int j=0;j<n_cosbin;j++){
      N_mom_cos_compare[j] +=N_mom_cos[i][j];
    }
  }

  TCanvas* c1 = new TCanvas("c1", "Lambda Polarization", 1000, 500);

  c1->Divide(2, 1);

  TMultiGraph *mg;
  TSpline5* spline[2][2];
  TGraphErrors *gr_exp[2];
  TGraphErrors *gr_model[2][2];
  TMultiGraph *mg_model[2];
  double cos_cm[n_cosbin];
  double cos_cm_err[n_cosbin];
  double model_value[2][2][n_cosbin];
  double model_value_err[2][2][n_cosbin];

  for (int mom = 0; mom < 2; mom++) {
    mg_model[mom] = new TMultiGraph();
    gr_exp[mom] = new TGraphErrors(n_cos,
                                    coscm,
                                    pol[mom],
                                    coscm_err,
                                    pol_err[mom]);
    c1->cd(mom + 1);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.15);
    c1->SetTopMargin(0.15);
    c1->SetBottomMargin(0.15);
    c1->SetGrid();
    gr_exp[mom]->SetMarkerStyle(20);
    gr_exp[mom]->GetXaxis()->SetTitle("cos#theta^{CM}_{#eta}");
    gr_exp[mom]->GetYaxis()->SetTitle("P_{#Lambda}");
    gr_exp[mom]->GetXaxis()->SetLimits(-1.0,1.0);
    gr_exp[mom]->GetYaxis()->SetRangeUser(-1.0, 1.0);
    gr_exp[mom]->Draw("AP");
    
      
    for (int theory = 0; theory < 2; ++theory) {
      TGraph* gr = new TGraph(cos_theta[mom][theory].size(),
			      &cos_theta[mom][theory][0],
			      &polarity[mom][theory][0]);
      gr->SetLineColor(theory == 0 ? colors[0] : colors[1]);
      
      spline[mom][theory] = new TSpline5(Form("spline_%d%d", mom,theory), gr);
      spline[mom][theory]->SetLineColor(gr->GetLineColor());
      //spline->SetLineStyle(2);
      spline[mom][theory]->SetLineWidth(5);
      spline[mom][theory]->Draw("same");
      for(int i=0;i<n_cosbin;i++){
	if(mom==0 && theory==0){
	  cos_cm[i] = xmin+step/2. + step*i;
	  cos_cm_err[i] = step/2.;
	}
	model_value[mom][theory][i] = spline[mom][theory]->Eval(cos_cm[i]);
	if(mom==0){
	  //model_value_err[mom][theory][i] = sqrt(3)/(alpha*sqrt((double)N_mom_cos[5][i]));
	  model_value_err[mom][theory][i] = sqrt(3)/(alpha*sqrt((double)N_mom_cos_compare[i]));
	}
	else if(mom==1){
	  model_value_err[mom][theory][i] = sqrt(3)/(alpha*sqrt((double)N_mom_cos[13][i]));
	}
	
      }
      gr_model[mom][theory] = new TGraphErrors(n_cosbin,cos_cm,model_value[mom][theory],cos_cm_err,model_value_err[mom][theory]);
      gr_model[mom][theory]->SetMarkerColor(theory == 0 ? colors[0] : colors[1]);
      gr_model[mom][theory]->SetLineColor(theory == 0 ? colors[0] : colors[1]);
      gr_model[mom][theory]->SetMarkerStyle(20);
      
      mg_model[mom]->Add(gr_model[mom][theory]);
      
    }
    mg_model[mom]->GetXaxis()->SetLimits(-1.0,1.0);
    mg_model[mom]->GetYaxis()->SetRangeUser(-1.0, 1.0);
    mg_model[mom]->GetXaxis()->SetTitle("cos#theta^{CM}_{#eta}");
    mg_model[mom]->GetYaxis()->SetTitle("P_{#Lambda}");
    mg_model[mom]->Add(gr_exp[mom]);
    gPad->SetGrid();
  }
  c1->cd(1);
  TLegend *le = new TLegend(0.6, 0.7, 0.88, 0.88);
  le->AddEntry(gr_exp[0],"Crystal Ball");
  le->AddEntry(spline[1][0],"D_{03}");
  le->AddEntry(spline[1][1],"P_{03}");
  le->Draw("same");


  TCanvas *c2= new TCanvas("c2");
  c2->Divide(2);

  for(int mom=0;mom<2;mom++){
    c2->cd(mom+1);
      gPad->SetGrid();
    mg_model[mom]->Draw("AP");

  }

  TGraphErrors *gr_polar_total = new TGraphErrors();
  for(int i=0;i<n_mombin;i++){
    //gr_polar_total->SetPoint(i,mom_kaon[i], sqrt(3)/(alpha*sqrt(N_mom_cos[i][7])));
    gr_polar_total->SetPoint(i,mom_kaon[i], sqrt(3)/(alpha*sqrt(N_mom[i])));
    gr_polar_total->SetPointError(i,mom_err,0);
    
  }

  TCanvas *c3 = new TCanvas("c3");
  c3->SetLeftMargin(0.15);
  c3->SetRightMargin(0.15);
  c3->SetTopMargin(0.15);
  c3->SetBottomMargin(0.15);
  gPad->SetGrid();
  gr_polar_total->SetMarkerStyle(20);
  //gr_polar_total->SetMarkerSize(1.5);
  gr_polar_total->SetTitle(";#it{p}^{#it{K}}_{Lab} [MeV/#it{c}];Statistical Uncertainty for P_{#Lambda}");
  gr_polar_total->Draw("AP");
  
  
  
  
}
