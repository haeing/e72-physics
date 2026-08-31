// lambda_p_2D.C
//
// Run:
//   root -l lambda_p_2D.C
//
// K- p -> eta Lambda @ pK = 735 MeV/c
// 2D circular LH2 target: diameter = 8 cm
//
// Estimate Lambda-p elastic scattering probability
// before Lambda decay / target exit.

#include <iostream>
#include <cmath>

#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLine.h"
#include "TEllipse.h"


// ============================================================
// Constants
// ============================================================

// masses [GeV/c^2]
const double mK      = 0.493677;
const double mp      = 0.938272;
const double meta    = 0.547862;
const double mLambda = 1.115683;

// beam momentum [GeV/c]
const double pK = 0.735;

// Lambda c*tau [cm]
const double ctauLambda = 7.845;


// ------------------------------------------------------------
// Target
// ------------------------------------------------------------

// diameter = 8 cm
const double targetRadius = 4.0; // cm


// ------------------------------------------------------------
// LH2
// ------------------------------------------------------------

const double rhoLH2 = 0.0708;       // g/cm^3
const double molarH2 = 2.01588;     // g/mol
const double NA = 6.02214076e23;

// number of protons / cm^3
const double nProton =
  rhoLH2 / molarH2 * NA * 2.0;


// ------------------------------------------------------------
// Lambda-p elastic cross section
//
// First approximation:
// sigma = 9 mb
// ------------------------------------------------------------

const double sigma_mb = 9.0;

// 1 mb = 1e-27 cm^2
const double sigma =
  sigma_mb * 1.0e-27;


// reconstructed Lambda -> p pi- events
const double nLambdaData = 1.7e6;


// ============================================================
// Distance to target boundary
//
// Target:
//
//       x^2 + z^2 = R^2
//
// Beam direction = +z
//
// Lambda starts at (x,z) and travels with direction
//
//       (ux,uz)
//
// ============================================================

double PathToBoundary(
		      double x,
		      double z,
		      double ux,
		      double uz)
{
  // Normalize just in case

  double norm =
    sqrt(ux*ux + uz*uz);

  ux /= norm;
  uz /= norm;


  // solve
  //
  // (x + L ux)^2 + (z + L uz)^2 = R^2
  //
  // L^2 + 2 b L + c = 0

  double b =
    x*ux + z*uz;

  double c =
    x*x + z*z
    - targetRadius*targetRadius;

  double D =
    b*b - c;

  if(D < 0)
    return 0.0;


  // positive solution

  double L =
    -b + sqrt(D);

  return L;
}


// ============================================================
// Main
// ============================================================

void yield()
{
  gStyle->SetOptStat(0);

  TRandom3 rnd(0);

  const long long NMC = 2000000;


  // ========================================================
  // Initial K- p system
  // ========================================================

  double EK =
    sqrt(pK*pK + mK*mK);


  TLorentzVector K(
		   0.0,
		   0.0,
		   pK,
		   EK
		   );


  TLorentzVector proton(
			0.0,
			0.0,
			0.0,
			mp
			);


  TLorentzVector initial =
    K + proton;


  double sqrtS =
    initial.M();


  TVector3 betaCM =
    initial.BoostVector();


  // ========================================================
  // Two-body CM momentum
  // ========================================================

  double s =
    sqrtS*sqrtS;


  double term1 =
    s - pow(mLambda + meta,2);

  double term2 =
    s - pow(mLambda - meta,2);


  double pStar =
    sqrt(term1*term2)
    / (2.0*sqrtS);


  double ELambdaStar =
    sqrt(
	 mLambda*mLambda
	 + pStar*pStar
	 );


  std::cout << std::endl;

  std::cout
    << "========================================"
    << std::endl;

  std::cout
    << "K- p -> eta Lambda"
    << std::endl;

  std::cout
    << "pK       = "
    << pK*1000
    << " MeV/c"
    << std::endl;

  std::cout
    << "sqrt(s)  = "
    << sqrtS*1000
    << " MeV"
    << std::endl;

  std::cout
    << "p*       = "
    << pStar*1000
    << " MeV/c"
    << std::endl;

  std::cout
    << "Target diameter = "
    << 2.0*targetRadius
    << " cm"
    << std::endl;

  std::cout
    << "n proton = "
    << nProton
    << " /cm3"
    << std::endl;

  std::cout
    << "========================================"
    << std::endl;


  // ========================================================
  // Histograms
  // ========================================================

  TH1D *hMomentum =
    new TH1D(
	     "hMomentum",
	     ";p_{#Lambda} [MeV/c];Counts",
	     200,300,700
	     );


  TH1D *hTheta =
    new TH1D(
	     "hTheta",
	     ";#theta_{#Lambda}^{lab} [deg];Counts",
	     180,0,180
	     );


  TH2D *hMomTheta =
    new TH2D(
	     "hMomTheta",
	     ";#theta_{#Lambda}^{lab} [deg];p_{#Lambda} [MeV/c]",
	     180,0,180,
	     200,300,700
	     );


  TH1D *hTargetPath =
    new TH1D(
	     "hTargetPath",
	     ";Path inside LH_{2} [cm];Counts",
	     200,0,8
	     );


  TH1D *hDecayLength =
    new TH1D(
	     "hDecayLength",
	     ";Mean #Lambda decay length [cm];Counts",
	     200,0,6
	     );


  TH1D *hEffectivePath =
    new TH1D(
	     "hEffectivePath",
	     ";Effective LH_{2} path [cm];Counts",
	     200,0,6
	     );


  TH1D *hScatterProb =
    new TH1D(
	     "hScatterProb",
	     ";P(#Lambda p #rightarrow #Lambda p) [%];Counts",
	     200,0,0.2
	     );


  TH2D *hVertex =
    new TH2D(
	     "hVertex",
	     ";x [cm];z [cm]",
	     100,-4,4,
	     100,-4,4
	     );


  // ========================================================
  // sums
  // ========================================================

  double sumPscatter = 0.0;
  double sumTargetPath = 0.0;
  double sumEffectivePath = 0.0;
  double sumDecayLength = 0.0;
  double sumExitProb = 0.0;


  // ========================================================
  // Event loop
  // ========================================================

  for(long long iev=0; iev<NMC; iev++)
    {

      // ====================================================
      // 1. Production point
      //
      // uniform in 2D circular target
      //
      // r = R sqrt(U)
      // ====================================================

      double r =
	targetRadius
	* sqrt(rnd.Uniform());


      double phiVertex =
	rnd.Uniform(0,2*TMath::Pi());


      double vx =
	r*cos(phiVertex);

      double vz =
	r*sin(phiVertex);


      hVertex->Fill(vx,vz);


      // ====================================================
      // 2. K- p -> eta Lambda
      //
      // isotropic Lambda in CM
      // ====================================================

      double cosThetaCM =
	rnd.Uniform(-1.0,1.0);


      double sinThetaCM =
	sqrt(
	     1.0
	     - cosThetaCM*cosThetaCM
	     );


      double phiCM =
	rnd.Uniform(
		    0,
		    2*TMath::Pi()
		    );


      double pxStar =
	pStar
	* sinThetaCM
	* cos(phiCM);


      double pyStar =
	pStar
	* sinThetaCM
	* sin(phiCM);


      double pzStar =
	pStar
	* cosThetaCM;


      TLorentzVector Lambda(
			    pxStar,
			    pyStar,
			    pzStar,
			    ELambdaStar
			    );


      // ====================================================
      // CM -> LAB
      // ====================================================

      Lambda.Boost(betaCM);


      double pLambda =
	Lambda.P();


      double thetaLab =
	Lambda.Theta();


      hMomentum->Fill(
		      pLambda*1000
		      );


      hTheta->Fill(
		   thetaLab
		   * TMath::RadToDeg()
		   );


      hMomTheta->Fill(
		      thetaLab*TMath::RadToDeg(),
		      pLambda*1000
		      );


      // ====================================================
      // 3. 2D Lambda direction
      //
      // Project Lambda onto x-z plane.
      // ====================================================

      double px =
	Lambda.Px();

      double pz =
	Lambda.Pz();


      double p2D =
	sqrt(
	     px*px
	     + pz*pz
	     );


      if(p2D < 1e-10)
	continue;


      double ux =
	px/p2D;

      double uz =
	pz/p2D;


      // ====================================================
      // 4. Path through target
      // ====================================================

      double Ltarget =
	PathToBoundary(
		       vx,
		       vz,
		       ux,
		       uz
		       );


      hTargetPath->Fill(
			Ltarget
			);


      // ====================================================
      // 5. Mean Lambda decay length
      //
      // L = beta gamma c tau
      //   = p/m c tau
      // ====================================================

      double Ldecay =
	pLambda
	/ mLambda
	* ctauLambda;


      hDecayLength->Fill(
			 Ldecay
			 );


      // ====================================================
      // 6. Effective path
      //
      // Account for Lambda decay:
      //
      // Leff =
      //
      // integral_0^Ltarget
      // exp(-L/Ldecay) dL
      //
      // =
      //
      // Ldecay *
      // (1-exp(-Ltarget/Ldecay))
      //
      // ====================================================

      double Leff =
	Ldecay
	* (
	   1.0
	   - exp(
		 -Ltarget/Ldecay
		 )
	   );


      hEffectivePath->Fill(
			   Leff
			   );


      // ====================================================
      // 7. Lambda-p elastic scattering probability
      // ====================================================

      double Pscatter =
	nProton
	* sigma
	* Leff;


      hScatterProb->Fill(
			 Pscatter*100.0
			 );


      // probability Lambda reaches target boundary
      // without decaying

      double Pexit =
	exp(
	    -Ltarget/Ldecay
            );


      // ====================================================
      // sums
      // ====================================================

      sumPscatter += Pscatter;

      sumTargetPath += Ltarget;

      sumEffectivePath += Leff;

      sumDecayLength += Ldecay;

      sumExitProb += Pexit;
    }


  // ========================================================
  // Results
  // ========================================================

  double meanPscatter =
    sumPscatter/NMC;


  double meanTargetPath =
    sumTargetPath/NMC;


  double meanEffectivePath =
    sumEffectivePath/NMC;


  double meanDecayLength =
    sumDecayLength/NMC;


  double meanExitProb =
    sumExitProb/NMC;


  double expectedEvents =
    nLambdaData
    * meanPscatter;


  std::cout << std::endl;

  std::cout
    << "=========== RESULTS ==========="
    << std::endl;


  std::cout
    << "Mean pLambda       = "
    << hMomentum->GetMean()
    << " MeV/c"
    << std::endl;


  std::cout
    << "Mean theta_lab     = "
    << hTheta->GetMean()
    << " deg"
    << std::endl;


  std::cout
    << "Mean target path   = "
    << meanTargetPath
    << " cm"
    << std::endl;


  std::cout
    << "Mean decay length  = "
    << meanDecayLength
    << " cm"
    << std::endl;


  std::cout
    << "Mean effective path = "
    << meanEffectivePath
    << " cm"
    << std::endl;


  std::cout
    << "Exit before decay  = "
    << meanExitProb*100
    << " %"
    << std::endl;


  std::cout
    << "Lambda-p elastic P = "
    << meanPscatter*100
    << " %"
    << std::endl;


  std::cout
    << std::endl;


  std::cout
    << "For "
    << nLambdaData
    << " reconstructed Lambdas:"
    << std::endl;


  std::cout
    << "Expected Lambda-p elastic events = "
    << expectedEvents
    << std::endl;


  std::cout
    << "==============================="
    << std::endl;


  // ========================================================
  // Draw
  // ========================================================

  TCanvas *c1 =
    new TCanvas(
		"c1",
		"Lambda kinematics",
		1000,800
		);


  c1->Divide(2,2);

  c1->cd(1);
  hMomentum->Draw();

  c1->cd(2);
  hTheta->Draw();

  c1->cd(3);
  hMomTheta->Draw("colz");

  c1->cd(4);
  hVertex->Draw("colz");


  // target circle

  TEllipse *circle =
    new TEllipse(
		 0,0,
		 targetRadius,
		 targetRadius
		 );

  circle->SetFillStyle(0);
  circle->SetLineWidth(2);
  circle->Draw("same");


  // --------------------------------------------------------

  TCanvas *c2 =
    new TCanvas(
		"c2",
		"Lambda path and scattering",
		1000,800
		);


  c2->Divide(2,2);

  c2->cd(1);
  hTargetPath->Draw();

  c2->cd(2);
  hDecayLength->Draw();

  c2->cd(3);
  hEffectivePath->Draw();

  c2->cd(4);
  hScatterProb->Draw();
}
