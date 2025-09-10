////////////////////////////////////////////////////////////////////////
//
// A class for dE/dx functions for likelihood based particleID
//
// sungbino@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef PHYSDEDX_H
#define PHYSDEDX_H

#include <map>
#include "Math/VavilovAccurate.h"
#include "TF1.h"

namespace pid {
  class PhysdEdx {

  public:

    PhysdEdx();
    PhysdEdx(int pdg);
    ~PhysdEdx();

    void SetPdgCode(int pdg);
    int GetPdgCode(){ return pdgcode;};

    double Landau_xi(double KE, double pitch);
    double Get_Wmax(double KE);
    double meandEdx(double KE);
    double MPVdEdx(double KE, double pitch);
    double KEtoMomentum(double KE);
    double MomentumtoKE(double momentum);

    bool dEdx_PDF(double KE, double pitch, double dEdx, double *PDF_y, double *PDF_maxy);
    double dEdx_Gaus_Sigma(double KE, double pitch);
  
  private:

    int pdgcode;
    double mass;
    int charge;

    double densityEffect(double beta, double gamma);

    double betaGamma(double KE);

    // == Bethe-Bloch parameters, https://indico.fnal.gov/event/14933/contributions/28526/attachments/17961/22583/Final_SIST_Paper.pdf
    const double rho = 1.39; // [g/cm3], density of LAr
    const double K = 0.307075; // [MeV cm2 / mol]
    const double Z = 18.; // atomic number of Ar
    const double A = 39.948; // [g / mol], atomic mass of Ar
    const double I = 197.0e-6; // [MeV], mean excitation energy, JINST 19 (2024) 01, P01009
    const double me = 0.511; // [Mev], mass of electron
    // == Parameters for the density correction
    const double density_C = 5.2146;
    const double density_y0 = 0.2;
    const double density_y1 = 3.0;
    const double density_a = 0.19559;
    const double density_k = 3.0;
  }; //
} // namespace 
#endif // PHYSDEDX_H
