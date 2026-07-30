void lambda_eta_phase() {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  gStyle->SetPalette(1);

  string outpdf = "lambda-eta-phase.pdf";
  
  Double_t Km = 0.493;
  Double_t pk = 0.735;
  Double_t Pm = 0.938;
  Double_t pipmm =0.1396;
  Double_t pi0m = 0.135;

  Double_t Lambdam = 1.115;
  Double_t etam = 0.548;
  TLorentzVector target(0.0, 0.0, 0.0, Pm);
  TLorentzVector beam(0.0, 0.0, pk, TMath::Sqrt(pk*pk+Km*Km));
  TLorentzVector W = beam + target;
 
   //(Momentum, Energy units are Gev/C, GeV)
  Double_t masses1[2] = {Lambdam, etam};
  Double_t masses2[2] = {Pm, pipmm};
  Double_t masses3[3] = { pipmm, pipmm, pi0m} ;
  Double_t masses4[3] = { pipmm, pipmm, 0.} ;
  
  
 
   TGenPhaseSpace event1;
   event1.SetDecay(W, 2, masses1);
 
   TH2F *h2 = new TH2F("h2","h2", 200,0.9,1, 200,0,1);
   TH2F *h3 = new TH2F("h3","h3", 200,0.9,1, 200,0,1);

   TH2F *h4 = new TH2F("h4","h4", 200,0.9,1, 200,0,1);
   TH2F *h5 = new TH2F("h5","h5", 200,0.,1, 200,0.,0.4);
   TH2F *h6 = new TH2F("h6","h6", 200,0.,1, 200,0.,0.4);
   TH2F *h7 = new TH2F("h7","h7", 200,0.,1, 200,0.,0.4);

   TH2F *h8 = new TH2F("h8","h8", 200,0.,1, 200,0.,0.7);

   const Int_t nParticles = 7;
   const Int_t nMomentumComponents = 4;
   const char *particleNames[nParticles] = {
      "lambda", "eta", "proton_from_lambda", "pion_from_lambda",
      "piplus_from_eta", "pi0_from_eta", "piminus_from_eta"
   };
   const char *particleTitles[nParticles] = {
      "#Lambda", "#eta", "p from #Lambda", "#pi from #Lambda",
      "#pi^{+} from #eta", "#pi^{0} from #eta", "#pi^{-} from #eta"
   };
   const char *componentNames[nMomentumComponents] = {"p", "px", "py", "pz"};
   const char *componentTitles[nMomentumComponents] = {
      "|#vec{p}|", "p_{x}", "p_{y}", "p_{z}"
   };
   TH1F *hMomentum[nParticles][nMomentumComponents];
   for (Int_t i = 0; i < nParticles; i++) {
      for (Int_t j = 0; j < nMomentumComponents; j++) {
         const Double_t xmin = (j == 0) ? 0.0 : -1.2;
         const Double_t xmax = 1.2;
         hMomentum[i][j] = new TH1F(
            Form("h_%s_%s", particleNames[i], componentNames[j]),
            Form("%s: %s in lab;%s [GeV/#it{c}];Events",
                 particleTitles[i], componentTitles[j], componentTitles[j]),
            240, xmin, xmax);
      }
   }


   TLorentzVector pLambda;
   TLorentzVector peta;


   TGenPhaseSpace event2;
   TGenPhaseSpace event3;
   TGenPhaseSpace event4;

 
   for (Int_t n=0;n<100000;n++) {
      Double_t weight1 = event1.Generate();

      TLorentzVector *pLambda_ori = event1.GetDecay(0);
      TLorentzVector *peta_ori = event1.GetDecay(1);

      pLambda.SetPxPyPzE(pLambda_ori->Px(),pLambda_ori->Py(),pLambda_ori->Pz(),pLambda_ori->E());
      peta.SetPxPyPzE(peta_ori->Px(),peta_ori->Py(),peta_ori->Pz(),peta_ori->E());
      

      event2.SetDecay(pLambda, 2, masses2);
      Double_t weight2 = event2.Generate();
      TLorentzVector *pP_lam = event2.GetDecay(0);
      TLorentzVector *pPi_lam = event2.GetDecay(1);
	
      h4->Fill(pP_lam->CosTheta(),pP_lam->P());

      h5->Fill(pPi_lam->CosTheta(),pPi_lam->P());
      
      event3.SetDecay(peta, 3, masses3);
      Double_t weight3 = event3.Generate();
	
      //      cout<<peta.P() << "\t" <<endl;
      TLorentzVector *pPip_eta = event3.GetDecay(0);
      //TLorentzVector *pPim_eta = event.GetDecay(1);
      TLorentzVector *pPi0_eta = event3.GetDecay(2);

      
      h6->Fill(pPip_eta->CosTheta(),pPip_eta->P());
      h7->Fill(pPi0_eta->CosTheta(),pPi0_eta->P());

      event4.SetDecay(peta, 3, masses4);
      Double_t weight4 = event4.Generate();

      TLorentzVector *pPim_eta2 = event4.GetDecay(0);
      h8->Fill(pPim_eta2->CosTheta(),pPim_eta2->P());

      TLorentzVector *labParticles[nParticles] = {
         pLambda_ori, peta_ori, pP_lam, pPi_lam,
         pPip_eta, pPi0_eta, pPim_eta2
      };
      for (Int_t i = 0; i < nParticles; i++) {
         hMomentum[i][0]->Fill(labParticles[i]->P());
         hMomentum[i][1]->Fill(labParticles[i]->Px());
         hMomentum[i][2]->Fill(labParticles[i]->Py());
         hMomentum[i][3]->Fill(labParticles[i]->Pz());
      }

      
      //delete pLambda;
      //delete peta;
	
      //delete event2;
      //delete event3;

      
 
 

      h2->Fill(pLambda_ori->CosTheta(),pLambda_ori->P());
      h3->Fill(peta_ori->CosTheta(),peta_ori->P());
 
   }
   TCanvas *c1 = new TCanvas("c1", "c1", 900, 700);
   TPaveText *p = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
   p->AddText("lambda-eta-phase.cc");
   TDatime now;
   p->AddText(Form("Generated at: %04d-%02d-%02d %02d:%02d:%02d",
		   now.GetYear(), now.GetMonth(), now.GetDay(),
		   now.GetHour(), now.GetMinute(), now.GetSecond()));
   p->Draw();
   c1->Print((outpdf + "(").c_str());
   c1->Clear();

  
   h2->SetTitle("#Lambda;cos#theta;momentum [GeV/#it{c}]");
   h2->Draw("colz");
   c1->Print(outpdf.c_str());
   c1->Clear();
   
   h3->SetMarkerColor(kBlack);
   h3->SetTitle("#eta;cos#theta;momentum [GeV/#it{c}]");
   h3->Draw("colz");
   c1->Print(outpdf.c_str());
   c1->Clear();
   
   h4->SetMarkerColor(kBlack);
   h4->SetTitle("P from #Lambda;cos#theta;momentum [GeV/#it{c}]");
   h4->Draw("colz");
   c1->Print(outpdf.c_str());
   c1->Clear();

   h5->SetMarkerColor(kBlack);
   h5->SetTitle("#pi from #Lambda;cos#theta;momentum [GeV/#it{c}]");
   h5->Draw("colz");
   c1->Print(outpdf.c_str());
   c1->Clear();

   h6->SetMarkerColor(kBlack);
   h6->SetTitle("#pi^+ from #eta;cos#theta;momentum [GeV/#it{c}]");
   h6->Draw("colz");
   c1->Print(outpdf.c_str());
   c1->Clear();

   h7->SetMarkerColor(kBlack);
   h7->SetTitle("#pi^0 from #eta;cos#theta;momentum [GeV/#it{c}]");
   h7->Draw("colz");
   c1->Print(outpdf.c_str());
   c1->Clear();

   h8->SetMarkerColor(kBlack);
   h8->SetTitle("#pi^- from #eta;cos#theta;momentum [GeV/#it{c}]");
   h8->Draw("colz");
   c1->Print(outpdf.c_str());
   c1->Clear();

   for (Int_t i = 0; i < nParticles; i++) {
      c1->Divide(2, 2);
      for (Int_t j = 0; j < nMomentumComponents; j++) {
         c1->cd(j + 1);
         hMomentum[i][j]->SetLineWidth(2);
         hMomentum[i][j]->Draw("hist");
      }
      if (i == nParticles - 1) {
         c1->Print((outpdf + ")").c_str());
      } else {
         c1->Print(outpdf.c_str());
      }
      c1->Clear();
   }
   
}
