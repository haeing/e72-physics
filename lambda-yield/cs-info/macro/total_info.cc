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

void convert_roots(double M1, double P1, double M2, double P2, double M1_err, double P1_err, double M2_err, double P2_err, double &roots, double &roots_err){
  double E1 = TMath::Sqrt(M1*M1 + P1*P1);
  double E2 = TMath::Sqrt(M2*M2 + P2*P2);
  double totalE = E1 + E2;
  double totalpz = P1 + P2;

  double error_E1 = TMath::Sqrt((M1 * M1 * M1_err * M1_err) / (E1 * E1) +
                                   (P1 * P1 * P1_err * P1_err) / (E1 * E1));
  double error_E2 = TMath::Sqrt((M2 * M2 * M2_err * M2_err) / (E2 * E2) +
                                   (P2 * P2 * P2_err * P2_err) / (E2 * E2));

  double error_totalE = TMath::Sqrt(error_E1 * error_E1 + error_E2 * error_E2);
  double error_totalpz = TMath::Sqrt(P1_err * P1_err + P2_err * P2_err);

  roots =  TMath::Sqrt(totalE*totalE - totalpz*totalpz);

  if(rootserr){
    roots_err =TMath::Sqrt((totalE * totalE * error_totalE * error_totalE) / (roots * roots) +
			 (totalpz * totalpz * error_totalpz * error_totalpz) / (roots * roots) );
  }
  else if(!rootserr){
    roots_err = 0;
  }


  /*
  double P1_low = P1 - P1_err/2;
  double P1_high = P1 + P1_err/2;

  double E1_low = TMath::Sqrt(M1*M1 + P1_low*P1_low);
  double E1_high = TMath::Sqrt(M1*M1 + P1_high*P1_high);

  double totalE_low = E1_low + E2;
  double totalE_high = E1_high + E2;
  double totalpz_low = P1_low + P2;
  double totalpz_high = P1_high + P2;
  

  double roots_low = TMath::Sqrt(totalE_low*totalE_low - totalpz_low*totalpz_low);
  double roots_high = TMath::Sqrt(totalE_high*totalE_high - totalpz_high*totalpz_high);

  roots_err = roots_high - roots_low;
  */
  
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


//Kp -> .. (DOI : 10.1016/0550-3213(70)90461-X)
static const string BC_doi = "10.1016/0550-3213(70)90461-X";
static const int BC_num = 19;
static const double BC_E[BC_num] = {1536, 1544, 1552, 1560, 1569, 1578, 1586, 1595, 1606, 1615, 1624, 1633, 1643, 1653, 1662, 1671, 1682, 1687, 1696};
static const double BC_E_err[BC_num] = {2.6, 2.8, 3.0, 3.1, 3.3, 3.5, 3.7, 3.8, 4.0, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.5, 3.6, 3.6, 3.7};

static const int BC_mode_num = 17;
static const string BC_mode[BC_mode_num] = {"#bar{K}^{0}n", "#bar{K}^{0}n#pi", "#bar{K}^{0}p#pi", "#Lambda#pi^{0}", "#Lambda#eta", "#Sigma^{0}#pi^{0}", "#Lambda+neutral", "#Lambda#pi^{+}#pi^{-}", "#Sigma^{0}#pi^{+}#pi^{-}", "#Lambda#pi^{+}#pi^{-}#pi^{0}", "#Sigma^{+}#pi^{-}", "#Sigma^{-}#pi^{+}", "#Sigma^{+}#pi^{-}#pi^{0}", "#Sigma^{-}#pi^{+}#pi^{0}", "K^{-}p", "K^{-}p#pi^{0}", "K^{-}#pi^{+}n"};
static const double BC_cs[BC_mode_num][BC_num]={
  {6.32, 5.26, 4.62, 4.42, 4.04, 4.28, 4.39, 3.99, 3.90, 3.55, 3.50, 3.85, 3.57, 2.88, 2.63, 2.45, 2.85, 3.41, 3.83},
  {-999, -999, -999, -999, -999, 0.00, 0.02, 0.04, 0.06, 0.06, 0.03, 0.01, 0.01, 0.06, 0.19, 0.12, 0.26, 0.15, 0.22},
  {-999, -999, -999, -999, -999, 0.00, 0.00, 0.00, 0.00, 0.02, 0.01, 0.00, 0.01, 0.04, 0.05, 0.07, 0.12, 0.07, 0.15},
  {3.52, 4.03, 2.97, 3.26, 2.69, 2.91, 3.09, 2.41, 2.57, 2.39, 2.86, 2.61, 2.75, 2.34, 3.08, 2.93, 2.90, 3.10, 3.14},
  {-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, 0.79, 0.80, 0.66, 0.59, 0.28},
  {4.43, 3.28, 2.72, 2.68, 1.86, 2.36, 1.74, 1.84, 2.07, 1.96, 1.48, 1.86, 1.88, 2.19, 2.29, 2.34, 1.98, 1.95, 1.52},
  {1.54, 0.78, 0.84, 0.76, 0.96, 0.78, 0.95, 1.03, 1.07, 0.89, 0.69, 0.97, 0.72, 0.55, 0.84, 0.82, 0.68, 1.10, 0.96},
  {1.48, 1.23, 1.35, 1.52, 1.62, 1.56, 1.77, 1.63, 1.82, 1.85, 1.73, 2.04, 2.12, 2.28, 2.80, 2.94, 3.02, 3.10, 3.44},
  {0.05, 0.02, 0.06, 0.06, 0.12, 0.12, 0.11, 0.16, 0.14, 0.20, 0.11, 0.21, 0.31, 0.31, 0.49, 0.63, 0.60, 0.62, 0.46},
  {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.01, 0.00, 0.01, 0.01, 0.09, 0.17, 0.23, 0.14, 0.09},
  {8.69, 6.34, 6.00, 5.18, 4.96, 5.82, 4.90, 4.29, 4.38, 3.92, 3.41, 3.60, 2.93, 3.06, 2.71, 2.45, 2.18, 2.60, 2.08},
  {4.62, 3.65, 3.45, 2.68, 2.47, 2.75, 2.79, 2.29, 2.33, 2.47, 2.39, 2.70, 3.07, 3.17, 4.01, 3.55, 3.57, 3.20, 2.32},
  {0.20, 0.12, 0.08, 0.16, 0.16, 0.36, 0.16, 0.17, 0.21, 0.16, 0.19, 0.27, 0.36, 0.42, 0.47, 0.72, 0.77, 0.82, 0.70},
  {0.06, 0.07, 0.11, 0.09, 0.08, 0.27, 0.24, 0.22, 0.24, 0.21, 0.27, 0.30, 0.29, 0.44, 0.51, 0.60, 0.67, 0.61, 0.63},
  {25.8, 23.8, 23.1, 21.5, 21.7, 19.3, 19.1, 17.5, 17.7, 18.6, 16.0, 16.4, 15.6, 15.2, 14.2, 15.9, 17.3, 18.6, 20.5},
  {-999, -999, -999, -999, -999, 0.00, 0.01, 0.00, 0.00, 0.02, 0.03, 0.03, 0.06, 0.02, 0.11, 0.13, 0.24, 0.30, 0.35},
  {-999, -999, -999, -999, -999, 0.00, 0.00, 0.02, 0.02, 0.00, 0.00 ,0.02, 0.00, 0.08, 0.11, 0.15, 0.09, 0.27, 0.27}
};

static const double BC_cs_err[BC_mode_num][BC_num]={
  {0.89, 0.47, 0.37, 0.40, 0.37, 0.44, 0.34, 0.34, 0.27, 0.27, 0.25, 0.27, 0.27, 0.24, 0.25, 0.21, 0.25, 0.27, 0.30},
  {-999, -999, -999, -999, -999, 0.02, 0.02, 0.04, 0.03, 0.03, 0.02, 0.01, 0.01, 0.03, 0.06, 0.04, 0.07, 0.05, 0.06},
  {-999, -999, -999, -999, -999, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.03, 0.03, 0.05, 0.03, 0.05},
  {0.70, 0.45, 0.30, 0.38, 0.37, 0.43, 0.34, 0.32, 0.25, 0.26, 0.24, 0.25, 0.26, 0.24, 0.31, 0.26, 0.30, 0.32, 0.33},
  {-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, 0.14, 0.08, 0.08, 0.07, 0.06},
  {0.68, 0.36, 0.27, 0.31, 0.24, 0.33, 0.23, 0.22, 0.19, 0.20, 0.16, 0.19, 0.19, 0.20, 0.24, 0.20, 0.21, 0.21, 0.20},
  {0.33, 0.21, 0.16, 0.19, 0.18, 0.22, 0.17, 0.18, 0.14, 0.14, 0.12, 0.13, 0.12, 0.12, 0.15, 0.13, 0.15, 0.17, 0.17},
  {0.26, 0.14, 0.13, 0.15, 0.16, 0.18, 0.15, 0.15, 0.13, 0.13, 0.12, 0.14, 0.15, 0.16, 0.20, 0.18, 0.20, 0.20, 0.23},
  {0.05, 0.02, 0.03, 0.03, 0.04, 0.04, 0.03, 0.04, 0.03, 0.04, 0.03, 0.04, 0.05, 0.05, 0.07, 0.07, 0.08, 0.08, 0.07},
  {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.03, 0.03, 0.04, 0.03, 0.03},
  {0.84, 0.41, 0.34, 0.34, 0.32, 0.41, 0.28, 0.27, 0.21, 0.20, 0.17, 0.18, 0.17, 0.18, 0.18, 0.14, 0.15, 0.17, 0.16},
  {0.49, 0.25, 0.20, 0.19, 0.18, 0.22, 0.18, 0.16, 0.13, 0.14, 0.13, 0.15, 0.17, 0.18, 0.24, 0.18, 0.21, 0.19, 0.15},
  {0.09, 0.04, 0.03, 0.04, 0.04, 0.07, 0.04, 0.04, 0.03, 0.03, 0.03, 0.04, 0.05, 0.05, 0.06, 0.07, 0.08, 0.08, 0.07},
  {0.04, 0.03, 0.03, 0.03, 0.03, 0.06, 0.04, 0.04, 0.04, 0.03, 0.04, 0.04, 0.04, 0.05, 0.06, 0.06, 0.07, 0.06, 0.07},
  {2.0,  1.7,  1.4,  1.3,  1.4,  1.3,  1.2,  1.1,  1.1,  1.1,  1.0,  1.0,  1.0,  1.0,  1.0,  1.1,  1.1,  1.6,  1.2 },
  {-999, -999, -999, -999, -999, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.03, 0.02, 0.06, 0.05, 0.08, 0.06, 0.07},
  {-999, -999, -999, -999, -999, 0.02, 0.01, 0.02, 0.02, 0.01, 0.01, 0.02, 0.02, 0.04, 0.06, 0.06, 0.05, 0.06, 0.06}
};

static const int color_map[20] = {1,2,3,4,5,6,7,9,11,12,15,20,28,30,34,38,40,46,51,65};
void total_info(){

  double Leta1_pK_err[Leta1_num];
  double Leta1_E[Leta1_num];
  double Leta1_E_err[Leta1_num];
  for(int i=0;i<Leta1_num;i++){
    Leta1_pK_err[i] = Leta1_pK_high[i] - Leta1_pK_low[i];
    convert_roots(mK, Leta1_pK[i], mp, 0, mK_err, Leta1_pK_err[i], mp_err, 0, Leta1_E[i], Leta1_E_err[i]);
  }

  double Leta2_E[Leta2_num];
  double Leta2_E_err[Leta2_num];
  for(int i=0;i<Leta2_num;i++){
    convert_roots(mK, Leta2_pK[i], mp, 0, mK_err, Leta2_pK_err[i], mp_err, 0, Leta2_E[i], Leta2_E_err[i]);
  }

  TGraphErrors *Leta1 = new TGraphErrors(Leta1_num, Leta1_E, Leta1_cs, Leta1_E_err, Leta1_cs_err);
  Leta1->SetMarkerStyle(21);
  Leta1->SetMarkerSize(0.5);
  Leta1->SetMarkerColor(color_map[1]);
  Leta1->SetLineColor(color_map[1]);
  Leta1->SetLineStyle(1);
  
  TGraphErrors *Leta2 = new TGraphErrors(Leta2_num, Leta2_E, Leta2_cs, Leta2_E_err, Leta2_cs_err);
  Leta2->SetMarkerStyle(22);
  Leta2->SetMarkerSize(0.5);
  Leta2->SetMarkerColor(color_map[1]);
  Leta2->SetLineColor(color_map[1]);
  Leta2->SetLineStyle(1);

  vector<vector<double>> BC_re_E;
  vector<vector<double>> BC_re_E_err;
  vector<vector<double>> BC_re_cs;
  vector<vector<double>> BC_re_cs_err;
  
  TGraphErrors *BC[BC_mode_num];
  for(int i=0;i<BC_mode_num;i++){
    BC_re_E.push_back(std::vector<double>());
    BC_re_E_err.push_back(std::vector<double>());
    BC_re_cs.push_back(std::vector<double>());
    BC_re_cs_err.push_back(std::vector<double>());
    
    for(int j=0;j<BC_num;j++){
      if(BC_cs[i][j]!=-999){
	BC_re_E[i].push_back(BC_E[j]);
	if(!rootserr)BC_re_E_err[i].push_back(0);
	else if(rootserr)BC_re_E_err[i].push_back(BC_E_err[j]);
	BC_re_cs[i].push_back(BC_cs[i][j]);
	BC_re_cs_err[i].push_back(BC_cs_err[i][j]);
      }
    }

    BC[i] = new TGraphErrors(BC_re_E[i].size(), BC_re_E[i].data(), BC_re_cs[i].data(), BC_re_E_err[i].data(), BC_re_cs_err[i].data());
    BC[i]->SetMarkerStyle(20);
    BC[i]->SetMarkerSize(0.5);
    BC[i]->SetMarkerColor(color_map[i+2]);
    BC[i]->SetLineColor(color_map[i+2]);
    BC[i]->SetLineStyle(1);
  }
  

  auto mg = new TMultiGraph();
  mg->Add(Leta1);
  mg->Add(Leta2);
  for(int i=0;i<BC_mode_num;i++){
    if(i!=6){
      mg->Add(BC[i]);
    }
  }

  mg->SetTitle("Cross Section Data;#sqrt{s} [MeV];#sigma [mb]");

  int x_min = 1530;
  int x_max = 1800;
  int y_min = 0;
  int y_max = 6;
  mg->GetXaxis()->SetRangeUser(x_min, x_max);
  mg->GetYaxis()->SetRangeUser(y_min, y_max);
  mg->Draw("ALP");

  TLegend *le = new TLegend(0.8,0.5,0.48,0.6);
  /*
  le->AddEntry(Leta1,Form("#Lambda#eta, %s",Leta1_doi.c_str()));
  le->AddEntry(Leta2,Form("#Lambda#eta, %s",Leta2_doi.c_str()));
  */
  le->AddEntry(Leta1,"#Lambda#eta");
  le->AddEntry(Leta2,"#Lambda#eta");
  for(int i=0;i<BC_mode_num;i++){
    if(i!=6)
      le->AddEntry(BC[i],Form("%s",BC_mode[i].c_str()));
  }
  

  le->Draw();

  //Mass Threshold
  TLine *threshold[BC_mode_num];
  double thre_m[BC_mode_num];
  thre_m[0] = mK0 + mn;
  thre_m[1] = mK0 + mn + mpi0;
  thre_m[2] = mK0 + mp + mpi;
  thre_m[3] = mLambda + mpi0;
  thre_m[4] = mLambda + meta;
  thre_m[5] = mSigma0 + mpi0;
  thre_m[6] = mLambda;
  thre_m[7] = mLambda + mpi + mpi;
  thre_m[8] = mSigma0 + mpi + mpi;

  auto axis = new TGaxis(x_min, y_max, x_max, y_max,x_min,x_max, 510, "-");
  //auto axis = mg->GetXaxis();
  for(int i=0;i<BC_mode_num;i++){
    threshold[i] = new TLine(thre_m[i], y_min, thre_m[i], y_max);
    threshold[i]->SetLineColor(color_map[i+2]);
    threshold[i]->SetLineStyle(2);
    threshold[i]->Draw();

    axis->ChangeLabelByValue(thre_m[i],-1,-1,-1,color_map[i+2],-1,Form("%s",BC_mode[i].c_str()));
    axis->Draw();
  }
  
  

  
}
