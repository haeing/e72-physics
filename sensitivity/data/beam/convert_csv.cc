void convert_csv(){

  std::string inputCsv   = "Nbeam_251123.csv";
  std::string outputRoot = "Nbeam_251123.root";

  std::ifstream fin(inputCsv);

  std::vector<double> p;
  std::vector<double> v;

  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;

    std::stringstream ss(line);
    std::string a, b;
    if (!std::getline(ss, a, ',')) continue;
    if (!std::getline(ss, b, ',')) continue;

    try {
      p.push_back(std::stod(a));
      v.push_back(std::stod(b));
    } catch (...) {
      std::cerr << "Warning: cannot parse line: " << line << std::endl;
      continue;
    }
  }

  int n = p.size();
  if (n < 2) {
    std::cerr << "Not enough data in CSV." << std::endl;
  }

 
  double step = p[1] - p[0];       // ex) 0.5
 
  std::vector<double> edges(n + 1);
  edges[0] = p[0] - 0.5 * step;
  for (int i = 1; i <= n; ++i) {
    edges[i] = edges[0] + i * step;
  }

 
  TFile *fout = new TFile(outputRoot.c_str(), "RECREATE");


  TH1D *hMom05 = new TH1D("hMom05",
			  "Momentum (0.5 step);Momentum;Value",
			  n, &edges[0]);

  for (int i = 0; i < n; ++i) {
    hMom05->SetBinContent(i + 1, v[i]);
  }

  int rebinFactor = static_cast<int>(std::round(2.0 / step));  // 2 / 0.5 = 4
  if (rebinFactor < 1) rebinFactor = 1;

  TH1D *hMom2 = (TH1D*) hMom05->Rebin(rebinFactor, "hMom2");

  hMom05->Write();
  hMom2->Write();

  fout->Close();


}
