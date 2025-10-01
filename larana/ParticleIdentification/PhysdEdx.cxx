////////////////////////////////////////////////////////////////////////
//
// A class for dE/dx functions for likelihood based particleID
//
// sungbino@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#include "PhysdEdx.h"
#include "TSpline.h"
#include <algorithm>
#include <cmath>
#include <iostream>

#include "cetlib/pow.h"

namespace {
  // == Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf
  const double rho = 1.39;   // [g/cm3], density of LAr
  const double K = 0.307075; // [MeV cm2 / mol]
  const double Z = 18.;      // atomic number of Ar
  const double A = 39.948;   // [g / mol], atomic mass of Ar
  const double I = 197.0e-6; // [MeV], mean excitation energy, JINST 19 (2024) 01, P01009
  const double me = 0.511;   // [Mev], mass of electron
  // == Parameters for the density correction
  const double density_C = 5.2146;
  const double density_y0 = 0.2;
  const double density_y1 = 3.0;
  const double density_a = 0.19559;
  const double density_k = 3.0;
}

pid::PhysdEdx::PhysdEdx(int pdg) : pdgcode(0), mass(0), charge(0)
{
  SetPdgCode(pdg);
}

void pid::PhysdEdx::SetPdgCode(int pdg)
{

  pdgcode = pdg;

  if (abs(pdgcode) == 13) { //muon
    mass = 105.6583755;
    charge = 1;
  }
  else if (abs(pdgcode) == 211) { //pion
    mass = 139.57039;
    charge = 1;
  }
  else if (abs(pdgcode) == 321) { //kaon
    mass = 493.677;
    charge = 1;
  }
  else if (pdgcode == 2212) { //proton
    mass = 938.27208816;
    charge = 1;
  }
  else {
    throw cet::exception("PhysdEdx") << "Unknown pdg code " << pdgcode;
  }
}

double pid::PhysdEdx::densityEffect(double beta, double gamma)
{
  // == Estimate the density correction
  double density_y = TMath::Log10(beta * gamma);
  double ln10 = TMath::Log(10);
  double this_delta = 0.;
  if (density_y > density_y1) { this_delta = 2.0 * ln10 * density_y - density_C; }
  else if (density_y < density_y0) {
    this_delta = 0.;
  }
  else {
    this_delta =
      2.0 * ln10 * density_y - density_C + density_a * pow(density_y1 - density_y, density_k);
  }

  return this_delta;
}

double pid::PhysdEdx::betaGamma(double KE)
{

  double gamma, beta;
  gamma = (KE + mass) / mass;
  beta = sqrt(1 - 1 / cet::square(gamma));

  return beta * gamma;
}

double pid::PhysdEdx::Landau_xi(double KE, double pitch)
{
  double gamma = (KE / mass) + 1.0;
  double beta = TMath::Sqrt(1 - (1.0 / (gamma * gamma)));
  double xi = rho * pitch * 0.5 * K * (Z / A) * cet::square(1. / beta);
  return xi;
}

double pid::PhysdEdx::Get_Wmax(double KE)
{
  double gamma = (KE / mass) + 1.0;
  double beta = TMath::Sqrt(1 - (1.0 / (gamma * gamma)));
  double Wmax = (2.0 * me * cet::square(beta * gamma)) /
                (1.0 + 2.0 * me * (gamma / mass) + cet::square(me / mass));

  return Wmax;
}

double pid::PhysdEdx::meandEdx(double KE)
{

  double gamma = (KE + mass) / mass;
  double beta = sqrt(1 - 1 / cet::square(gamma));
  double wmax = Get_Wmax(KE);
  double dEdX =
    (rho * K * Z * cet::square(charge)) / (A * cet::square(beta)) *
    (0.5 * log(2 * me * cet::square(gamma) * cet::square(beta) * wmax / cet::square(I)) -
     cet::square(beta) - densityEffect(beta, gamma) / 2);

  return dEdX;
}

double pid::PhysdEdx::MPVdEdx(double KE, double pitch)
{

  //KE is kinetic energy in MeV
  //pitch is in cm
  double gamma = (KE + mass) / mass;
  double beta = sqrt(1 - 1 / cet::square(gamma));

  double xi = Landau_xi(KE, pitch);

  double eloss_mpv = xi *
                     (log(2 * me * cet::square(gamma) * cet::square(beta) / I) + log(xi / I) + 0.2 -
                      cet::square(beta) - densityEffect(beta, gamma)) /
                     pitch;

  return eloss_mpv;
}

double pid::PhysdEdx::KEtoMomentum(double KE)
{
  return sqrt(cet::square(KE) + 2.0 * KE * mass);
}

double pid::PhysdEdx::MomentumtoKE(double momentum)
{
  return sqrt(cet::square(momentum) + cet::square(mass)) - mass;
}

bool pid::PhysdEdx::dEdx_PDF(double KE, double pitch, double dEdx, double* PDF_y, double* PDF_maxy)
{
  auto pdf_lambda = [this](double* x, double* par) -> double {
    double a = par[2] / par[4];
    double b = (0.422784 + par[1] + log(par[0])) * par[2] / par[4] + par[3];
    double y = (x[0] - b) / a;

    double this_vav = 0.;

    if (par[0] < 0.01) { // == Landau
      this_vav = TMath::Landau(y);
      this_vav = this_vav / a;
    }
    else if (par[0] > 10.) { // == Gaussian
      double mu = vav.Mean(par[0], par[1]);
      double sigma = sqrt(vav.Variance(par[0], par[1]));
      this_vav = TMath::Gaus(y, mu, sigma);
    }
    else { // == Vavilov
      this_vav = vav.Pdf(y, par[0], par[1]);
      this_vav = this_vav / a;
    }

    return this_vav;
  };

  double gamma = (KE / mass) + 1.0;
  double beta = TMath::Sqrt(1 - (1.0 / (gamma * gamma)));
  double this_xi = Landau_xi(KE, pitch);
  double this_Wmax = Get_Wmax(KE);
  double this_kappa = this_xi / this_Wmax;
  double this_dEdx_BB = meandEdx(KE);
  double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch};

  TF1 PDF("", pdf_lambda, 0., 50., 5);
  PDF.SetParameters(par[0], par[1], par[2], par[3], par[4]);
  double this_PDF_y = PDF.Eval(dEdx);
  PDF_y[0] = this_PDF_y;

  if (par[0] > 0.01 && par[0] < 10.) {
    double mu = vav.Mean(par[0], par[1]);
    double sigma = sqrt(vav.Variance(par[0], par[1]));
    TF1 PDF_narrow("PDF_narrow", pdf_lambda, mu - sigma, mu + sigma, 5);
    PDF_narrow.SetParameters(par[0], par[1], par[2], par[3], par[4]);
    PDF_maxy[0] = PDF_narrow.GetMaximum();
  }
  else {
    PDF_maxy[0] = PDF.GetMaximum();
  }

  return this_PDF_y > 0.;
}

double pid::PhysdEdx::dEdx_Gaus_Sigma(double KE, double pitch)
{

  double gamma = (KE / mass) + 1.0;
  double beta = TMath::Sqrt(1 - (1.0 / (gamma * gamma)));
  double this_xi = Landau_xi(KE, pitch);
  double this_Wmax = Get_Wmax(KE);
  double this_kappa = this_xi / this_Wmax;

  double sigma = sqrt(vav.Variance(this_kappa, beta * beta));

  return sigma;
}
