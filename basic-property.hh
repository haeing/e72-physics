#ifndef BASIC_PROPERTY_HH
#define BASIC_PROPERTY_HH

#include "TMath.h"

// ============================================================
// Particle masses [MeV/c^2]
// ============================================================
static const double mK       = 493.677;
static const double mp       = 938.27208816;
static const double mn       = 939.5654205;
static const double meta     = 547.862;
static const double mLambda  = 1115.683;
static const double mK0      = 497.611;
static const double mSigmap  = 1189.37;
static const double mSigma0  = 1192.642;
static const double mSigmam  = 1197.449;
static const double mpi      = 139.57039;
static const double mpi0     = 134.9768;

// ============================================================
// Particle properties
// ============================================================
const double tauK = 1.238e-8;      // s

// ============================================================
// Beam intensity / running condition
// ============================================================
static const double spill_length = 4.24; // s

// distance between BAC and target center [m]
static const double dist_bac_tgt = 0.4983;


static const double pK_intensity[3] = {700, 735, 750};       // MeV/c
static const double intensity[3]    = {27400, 38400, 44000}; // /spill

static const double FK = 3.8e4; // beam particles per spill

static const double T_scan = 0.5 * 24.0 * 3600.0 / spill_length;
static const double T_fix  = 5.5 * 24.0 * 3600.0 / spill_length;

static const double eff_acc = 1.0;

// representative beam momentum [MeV/c]
static const double pK_fix = 735.0;

// K- lifetime [s]
static const double K_tau = 1.238e-8;

// beam-profile ratio from target upstream face to cross-sectional area
static const double Rbeam = 0.69;

// ============================================================
// Useful conversion functions
// target p at rest, K- beam in lab
// ============================================================
inline double SqrtSFromPK(double pK)
{
  double EK = TMath::Sqrt(mK*mK + pK*pK);
  double s  = mK*mK + mp*mp + 2.0 * mp * EK;
  return TMath::Sqrt(s);
}

inline double PKFromSqrtS(double sqrts)
{
  double s = sqrts * sqrts;
  double EK = (s - mK*mK - mp*mp) / (2.0 * mp);
  return TMath::Sqrt(EK*EK - mK*mK);
}

inline double PKToSqrtSErr(double pK, double pKerr)
{
  double e1 = SqrtSFromPK(pK + pKerr);
  double e2 = SqrtSFromPK(pK - pKerr);
  return 0.5 * TMath::Abs(e1 - e2);
}

inline double SqrtSToPKErr(double sqrts, double sqrts_err)
{
  double p1 = PKFromSqrtS(sqrts + sqrts_err);
  double p2 = PKFromSqrtS(sqrts - sqrts_err);
  return 0.5 * TMath::Abs(p1 - p2);
}

inline double KaonSurvivalProb(double pK_MeV, double L_m)
{
  const double c = 299792458.0;      // m/s

  double beta_gamma = pK_MeV / mK;
  double decay_length = beta_gamma * c * tauK; // m

  return TMath::Exp(-L_m / decay_length);
}

#endif
