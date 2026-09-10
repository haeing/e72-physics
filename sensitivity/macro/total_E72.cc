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
static const double mpi0 = 134.9768;       //(uubar - ddbar)/sqrt(2)


static const bool rootserr = false;

double pK_to_roots(double *x,double *p){
  double pK = x[0];
  return sqrt(pow(sqrt(mK*mK+pK*pK)+mp,2)-pK*pK);
}


//Energy : MeV
//Cross section : mb


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



void total_E72(){

  double Leta1_pK_err[Leta1_num];
  //double Leta1_E[Leta1_num];
  //double Leta1_E_err[Leta1_num];
  for(int i=0;i<Leta1_num;i++){
    Leta1_pK_err[i] = Leta1_pK_high[i] - Leta1_pK_low[i];
    //convert_roots(mK, Leta1_pK[i], mp, 0, mK_err, Leta1_pK_err[i], mp_err, 0, Leta1_E[i], Leta1_E_err[i]);
  }
  /*
  double echeck;
  double ee;
  convert_roots(mK,670,mp,0,mK_err,0,mp_err,0,echeck,ee);

  double Leta2_E[Leta2_num];
  double Leta2_E_err[Leta2_num];
  for(int i=0;i<Leta2_num;i++){
    convert_roots(mK, Leta2_pK[i], mp, 0, mK_err, Leta2_pK_err[i], mp_err, 0, Leta2_E[i], Leta2_E_err[i]);
  }
  */

  //TGraphErrors *Leta1 = new TGraphErrors(Leta1_num, Leta1_E, Leta1_cs, Leta1_E_err, Leta1_cs_err);
  TGraphErrors *Leta1 = new TGraphErrors(Leta1_num, Leta1_pK, Leta1_cs, Leta1_pK_err, Leta1_cs_err);
  Leta1->SetMarkerStyle(21);
  Leta1->SetMarkerSize(1);
  Leta1->SetLineStyle(1);
  
  //TGraphErrors *Leta2 = new TGraphErrors(Leta2_num, Leta2_E, Leta2_cs, Leta2_E_err, Leta2_cs_err);
  TGraphErrors *Leta2 = new TGraphErrors(Leta2_num, Leta2_pK, Leta2_cs, Leta2_pK_err, Leta2_cs_err);
  Leta2->SetMarkerStyle(22);
  Leta2->SetMarkerSize(1);
  Leta2->SetLineStyle(1);

  TCanvas* c1 = new TCanvas("c1","c1");
  c1->cd();
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  c1->SetTopMargin(0.15);
  c1->SetBottomMargin(0.15);

  auto mg = new TMultiGraph();
  mg->Add(Leta1);
  mg->Add(Leta2);
  
  mg->SetTitle("Cross Section Data;#it{p_{K}} [MeV/#it{c}];#sigma [mb]");
  
  int x_min = 710;
  int x_max = 800;
  int y_min = 0;
  int y_max = 2;
  mg->GetXaxis()->SetLimits(x_min, x_max);
  mg->GetYaxis()->SetRangeUser(y_min, y_max);
  mg->Draw("AP");

  TLegend *le = new TLegend(0.8,0.5,0.48,0.6);

  le->AddEntry(Leta1,"A. Starostin #it{et al.}");
  le->AddEntry(Leta2,"R. Rader #it{et al.}");
  le->Draw();
  
  //TF1* f1 = new TF1("f1", pK_to_roots, x_min, x_max, 0);
  cout<<sqrt(pow(sqrt(pow(493.677,2)+pow(735,2))+938.272,2)-735*735)<<endl;
  TF1* f1 = new TF1("f1", "sqrt((sqrt(493.677^2 + x^2) + 938.272)^2 - x^2)", x_min, x_max);
  TGaxis* axis = new TGaxis(x_min, y_max, x_max, y_max, "f1", 510, "-", 0.03);
  
  //axis->Draw("same");

}
