#ifndef TPC_EVENT_DISPLAY_HELPER_HH
#define TPC_EVENT_DISPLAY_HELPER_HH

#include <TH2Poly.h>
#include <TMath.h>

#include "TPCPadHelper.hh"

namespace tpcdisp
{
  inline TH2Poly* MakeTPCPadMap(const char* name,
                                const char* title = "",
                                double zOffset = ZTarget)
  {
    auto h = new TH2Poly(name,title,MinZ,MaxZ,MinX,MaxX);
    h->SetName(name);
    h->SetTitle(title);
    h->GetXaxis()->SetTitle("Z [mm]");
    h->GetYaxis()->SetTitle("X [mm]");
    Double_t X[5];
    Double_t Y[5];

    for (Int_t l = 0; l < NumOfLayersTPC; ++l) {
      Double_t pLength = tpc::padParameter[l][5];
      Double_t st =
        (180. - (360. / tpc::padParameter[l][3]) * tpc::padParameter[l][1] / 2.);

      Double_t sTheta = (-1. + st / 180.) * TMath::Pi();
      Double_t dTheta = (360. / tpc::padParameter[l][3]) / 180. * TMath::Pi();
      Double_t cRad   = tpc::padParameter[l][2];
      Int_t    nPad   = tpc::padParameter[l][1];

      for (Int_t j = 0; j < nPad; ++j) {
        X[1] = (cRad + pLength / 2.) * TMath::Cos(j * dTheta + sTheta);
        X[2] = (cRad + pLength / 2.) * TMath::Cos((j + 1) * dTheta + sTheta);
        X[3] = (cRad - pLength / 2.) * TMath::Cos((j + 1) * dTheta + sTheta);
        X[4] = (cRad - pLength / 2.) * TMath::Cos(j * dTheta + sTheta);
        X[0] = X[4];

        Y[1] = (cRad + pLength / 2.) * TMath::Sin(j * dTheta + sTheta);
        Y[2] = (cRad + pLength / 2.) * TMath::Sin((j + 1) * dTheta + sTheta);
        Y[3] = (cRad - pLength / 2.) * TMath::Sin((j + 1) * dTheta + sTheta);
        Y[4] = (cRad - pLength / 2.) * TMath::Sin(j * dTheta + sTheta);
        Y[0] = Y[4];

        for (Int_t k = 0; k < 5; ++k) X[k] += zOffset;

        h->AddBin(5, X, Y);
      }
    }

    return h;
  }

  inline TEllipse* DrawTarget(double zTarget = ZTarget,
			      double radius  = 40.0,
			      Color_t color  = kRed,
			      Style_t style  = 1,
			      Width_t width  = 2)
  {
    auto target = new TEllipse(zTarget, 0.0, radius, radius);
    target->SetFillStyle(0);
    target->SetLineColor(color);
    target->SetLineStyle(style);
    target->SetLineWidth(width);
    target->Draw("same");
    return target;
  }

  inline TEllipse* DrawTargetHolder(double zTarget = ZTarget,
				    double radius  = 50.0,
				    Color_t color  = kGreen,
				    Style_t style  = 1,
				    Width_t width  = 2)
  {
    auto target = new TEllipse(zTarget, 0.0, radius, radius);
    target->SetFillStyle(0);
    target->SetLineColor(color);
    target->SetLineStyle(style);
    target->SetLineWidth(width);
    target->Draw("same");
    return target;
  }
  
}

#endif
