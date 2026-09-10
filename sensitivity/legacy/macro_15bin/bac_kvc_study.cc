static const double spill_length = 4.24; // [s]

static const double pedestal_hz = 2000; //Hz
void bac_kvc_study(){
  

  double bac_intensity[4] = {165.938, 546.440, 694.715, 913.692};


   
  double bac_eff[4] = {0.998, 0.998, 0.998, 0.998};

  double kvc_HV[5] = {58, 59, 60, 61, 62};
  double kvc_ineff[5] = {0.017, 0.074, 0.300, 0.681, 0.948};
  double kvc_noise[5];
  for(int i=0;i<5;i++)kvc_noise[i] = kvc_ineff[i] * pedestal_hz;
  
  double kvc_eff[5] = {0.487, 0.652, 0.808, 0.938, 0.991};

  

  TGraph *bac = new TGraph(4, bac_intensity, bac_eff);
  TGraph *g_kvc_eff = new TGraph(5, kvc_HV, kvc_eff);
  TGraph *g_kvc_ineff = new TGraph(5, kvc_HV, kvc_ineff);
  TGraph *g_kvc_noise = new TGraph(5, kvc_HV, kvc_noise);

  TMultiGraph *g_kvc = new TMultiGraph();
  g_kvc->Add(g_kvc_eff);
  g_kvc->Add(g_kvc_ineff);
  //g_kvc->Add(g_kvc_noise);

  TCanvas *c1 = new TCanvas("c1","c1");
  bac->SetMarkerStyle(20);
  bac->SetTitle(";Intensity [k / spill]; Efficiency");
  bac->GetYaxis()->SetRangeUser(0.5,1.05);
  bac->Draw("AP");

  TCanvas *c2 = new TCanvas("c2","c2");
  //c2->SetLogy();
  g_kvc_eff->SetMarkerStyle(20);
  g_kvc_ineff->SetMarkerStyle(20);
  g_kvc_ineff->SetMarkerColor(kRed);

  g_kvc_noise->SetMarkerStyle(20);
  g_kvc_noise->SetMarkerColor(kRed);

  double y_min = 0.01;  // 왼쪽 Y축 최소값
  double y_max = 1.05 ; // 왼쪽 Y축 최대값
  g_kvc->SetTitle(";MPPC HV [V]; Efficiency");
  g_kvc->GetXaxis()->SetLimits(57.5,62.5);
  g_kvc->GetYaxis()->SetRangeUser(y_min,y_max);

  g_kvc_eff->SetTitle(";MPPC HV [V]; Efficiency");
  g_kvc_eff->GetXaxis()->SetLimits(57.5,62.5);
  
  g_kvc_eff->Draw("AP");

  


  TGaxis* rightAxis = new TGaxis(62.5,0,62.5,1.05,0,1.05,510, "+L");
  rightAxis->SetLineColor(kRed);
  rightAxis->SetLabelColor(kRed);
  rightAxis->SetTitle("Inefficiency");
  rightAxis->SetTitleColor(kRed);
  //rightAxis->Draw();

  TCanvas *c3 = new TCanvas("c3","c3");
  c3->SetLogy();
  g_kvc_noise->SetTitle(";MPPC HV [V]; Noise Rate [Hz]");
  g_kvc_noise->GetXaxis()->SetLimits(57.5,62.5);
  g_kvc_noise->GetYaxis()->SetRangeUser(1,5000);
  g_kvc_noise->Draw("AP");

  
  
}
