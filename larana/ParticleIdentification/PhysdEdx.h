////////////////////////////////////////////////////////////////////////
//
// A class for dE/dx functions for likelihood based particleID
//
// sungbino@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef PHYSDEDX_H
#define PHYSDEDX_H

#include "Math/VavilovAccurate.h"
#include "TF1.h"
#include <map>

namespace pid {
  class PhysdEdx {

  public:
    PhysdEdx(int pdg);
    int GetPdgCode() { return pdgcode; };

    double Landau_xi(double KE, double pitch);
    double Get_Wmax(double KE);
    double meandEdx(double KE);
    double MPVdEdx(double KE, double pitch);
    double KEtoMomentum(double KE);
    double MomentumtoKE(double momentum);

    bool dEdx_PDF(double KE, double pitch, double dEdx, double* PDF_y, double* PDF_maxy);
    double dEdx_Gaus_Sigma(double KE, double pitch);

  private:
    ROOT::Math::VavilovAccurate vav;
    void SetPdgCode(int pdg);
    int pdgcode;
    double mass;
    int charge;

    double densityEffect(double beta, double gamma);
    double betaGamma(double KE);
  }; //
} // namespace
#endif // PHYSDEDX_H
