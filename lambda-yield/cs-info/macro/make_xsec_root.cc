#include "TFile.h"
#include "TGraphErrors.h"
#include "TNamed.h"
#include "TString.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

using namespace std;

#include "total_info.cc"

struct Point {
  double x;
  double y;
  double ex;
  double ey;
};

double PkFromSqrtS(double sqrts)
{
  double s = sqrts * sqrts;

  double E_K_lab = (s - mK*mK - mp*mp) / (2.0 * mp);

  double pK = TMath::Sqrt(E_K_lab*E_K_lab - mK*mK);

  return pK;
}

double SqrtSFromPk(double pK)
{
  double EK = TMath::Sqrt(mK*mK + pK*pK);
  double s = mK*mK + mp*mp + 2.0 * mp * EK;

  return TMath::Sqrt(s);
}

double ConvertPkErrToSqrtSErr(double pK, double pKerr)
{
  double s0 = SqrtSFromPk(pK);
  double s1 = SqrtSFromPk(pK + pKerr);
  double s2 = SqrtSFromPk(pK - pKerr);

  return 0.5 * fabs(s1 - s2);
}

double ConvertSqrtSErrToPkErr(double sqrts, double sqrts_err)
{
  double p1 = PkFromSqrtS(sqrts + sqrts_err);
  double p2 = PkFromSqrtS(sqrts - sqrts_err);

  return 0.5 * fabs(p1 - p2);
}

void AddInfo(TGraphErrors *g,
             const char *reaction,
             const char *doi,
             const char *x_type)
{
  g->GetListOfFunctions()->Add(new TNamed("Reaction", reaction));
  g->GetListOfFunctions()->Add(new TNamed("DOI", doi));
  g->GetListOfFunctions()->Add(new TNamed("XType", x_type));
}

void SaveGraph(TFile *f,
               const char *name,
               const char *title,
               vector<Point> pts,
               const char *reaction,
               const char *doi,
               const char *x_title,
               const char *x_type)
{
  sort(pts.begin(), pts.end(),
       [](const Point &a, const Point &b){ return a.x < b.x; });

  int n = pts.size();

  vector<double> x(n), y(n), ex(n), ey(n);

  for (int i = 0; i < n; i++) {
    x[i]  = pts[i].x;
    y[i]  = pts[i].y;
    ex[i] = pts[i].ex;
    ey[i] = pts[i].ey;
  }

  TGraphErrors *g =
    new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());

  g->SetName(name);
  g->SetTitle(Form("%s;%s;#sigma [mb]", title, x_title));

  AddInfo(g, reaction, doi, x_type);

  f->cd();
  g->Write();
}

void make_xsec_root()
{
  TFile *fout = new TFile("total_cross_section.root", "RECREATE");

  // ============================================================
  // Lambda eta : Leta1 + Leta2 combined
  // ============================================================
  vector<Point> Lambda_eta_E;
  vector<Point> Lambda_eta_pK;

  for (int i = 0; i < Leta1_num; i++) {

    double pKerr =
      0.5 * (Leta1_pK_high[i] - Leta1_pK_low[i]);

    double E = SqrtSFromPk(Leta1_pK[i]);
    double Eerr = ConvertPkErrToSqrtSErr(Leta1_pK[i], pKerr);

    Lambda_eta_pK.push_back({
      Leta1_pK[i],
      Leta1_cs[i],
      pKerr,
      Leta1_cs_err[i]
    });

    Lambda_eta_E.push_back({
      E,
      Leta1_cs[i],
      Eerr,
      Leta1_cs_err[i]
    });
  }

  for (int i = 0; i < Leta2_num; i++) {

    double E = SqrtSFromPk(Leta2_pK[i]);
    double Eerr = ConvertPkErrToSqrtSErr(Leta2_pK[i], Leta2_pK_err[i]);

    Lambda_eta_pK.push_back({
      Leta2_pK[i],
      Leta2_cs[i],
      Leta2_pK_err[i],
      Leta2_cs_err[i]
    });

    Lambda_eta_E.push_back({
      E,
      Leta2_cs[i],
      Eerr,
      Leta2_cs_err[i]
    });
  }

  SaveGraph(fout,
            "g_Lambda_eta_E",
            "K^{-}p #rightarrow #Lambda#eta",
            Lambda_eta_E,
            "K- p -> Lambda eta",
            Form("%s, %s", Leta1_doi.c_str(), Leta2_doi.c_str()),
            "#sqrt{s} [MeV]",
            "sqrt_s");

  SaveGraph(fout,
            "g_Lambda_eta_pK",
            "K^{-}p #rightarrow #Lambda#eta",
            Lambda_eta_pK,
            "K- p -> Lambda eta",
            Form("%s, %s", Leta1_doi.c_str(), Leta2_doi.c_str()),
            "p_{K^{-}} [MeV/c]",
            "pK_lab");

  // ============================================================
  // Bubble chamber data
  // ============================================================
  const char *graph_base[BC_mode_num] = {
    "g_K0bar_n",
    "g_K0bar_n_pi",
    "g_K0bar_p_pi",
    "g_Lambda_pi0",
    "g_Lambda_eta_BC",
    "g_Sigma0_pi0",
    "g_Lambda_neutral",
    "g_Lambda_pip_pim",
    "g_Sigma0_pip_pim",
    "g_Lambda_pip_pim_pi0",
    "g_Sigmap_pim",
    "g_Sigmam_pip",
    "g_Sigmap_pim_pi0",
    "g_Sigmam_pip_pi0",
    "g_Km_p",
    "g_Km_p_pi0",
    "g_Km_pip_n"
  };

  const char *title[BC_mode_num] = {
    "K^{-}p #rightarrow #bar{K}^{0}n",
    "K^{-}p #rightarrow #bar{K}^{0}n#pi",
    "K^{-}p #rightarrow #bar{K}^{0}p#pi",
    "K^{-}p #rightarrow #Lambda#pi^{0}",
    "K^{-}p #rightarrow #Lambda#eta",
    "K^{-}p #rightarrow #Sigma^{0}#pi^{0}",
    "K^{-}p #rightarrow #Lambda + neutral",
    "K^{-}p #rightarrow #Lambda#pi^{+}#pi^{-}",
    "K^{-}p #rightarrow #Sigma^{0}#pi^{+}#pi^{-}",
    "K^{-}p #rightarrow #Lambda#pi^{+}#pi^{-}#pi^{0}",
    "K^{-}p #rightarrow #Sigma^{+}#pi^{-}",
    "K^{-}p #rightarrow #Sigma^{-}#pi^{+}",
    "K^{-}p #rightarrow #Sigma^{+}#pi^{-}#pi^{0}",
    "K^{-}p #rightarrow #Sigma^{-}#pi^{+}#pi^{0}",
    "K^{-}p #rightarrow K^{-}p",
    "K^{-}p #rightarrow K^{-}p#pi^{0}",
    "K^{-}p #rightarrow K^{-}#pi^{+}n"
  };

  const char *reaction[BC_mode_num] = {
    "K- p -> K0bar n",
    "K- p -> K0bar n pi",
    "K- p -> K0bar p pi",
    "K- p -> Lambda pi0",
    "K- p -> Lambda eta",
    "K- p -> Sigma0 pi0",
    "K- p -> Lambda neutral",
    "K- p -> Lambda pi+ pi-",
    "K- p -> Sigma0 pi+ pi-",
    "K- p -> Lambda pi+ pi- pi0",
    "K- p -> Sigma+ pi-",
    "K- p -> Sigma- pi+",
    "K- p -> Sigma+ pi- pi0",
    "K- p -> Sigma- pi+ pi0",
    "K- p -> K- p",
    "K- p -> K- p pi0",
    "K- p -> K- pi+ n"
  };

  for (int i = 0; i < BC_mode_num; i++) {

    vector<Point> pts_E;
    vector<Point> pts_pK;

    for (int j = 0; j < BC_num; j++) {

      if (BC_cs[i][j] == -999) continue;

      double pK = PkFromSqrtS(BC_E[j]);
      double pKerr = ConvertSqrtSErrToPkErr(BC_E[j], BC_E_err[j]);

      pts_E.push_back({
        BC_E[j],
        BC_cs[i][j],
        BC_E_err[j],
        BC_cs_err[i][j]
      });

      pts_pK.push_back({
        pK,
        BC_cs[i][j],
        pKerr,
        BC_cs_err[i][j]
      });
    }

    if (pts_E.empty()) continue;

    SaveGraph(fout,
              Form("%s_E", graph_base[i]),
              title[i],
              pts_E,
              reaction[i],
              BC_doi.c_str(),
              "#sqrt{s} [MeV]",
              "sqrt_s");

    SaveGraph(fout,
              Form("%s_pK", graph_base[i]),
              title[i],
              pts_pK,
              reaction[i],
              BC_doi.c_str(),
              "p_{K^{-}} [MeV/c]",
              "pK_lab");
  }

  fout->Close();

  cout << "Saved: total_cross_section.root" << endl;
}
