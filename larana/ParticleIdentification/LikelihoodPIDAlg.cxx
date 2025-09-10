////////////////////////////////////////////////////////////////////////
//
// A likelihood based particleID
//
// sungbino@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#include "larana/ParticleIdentification/LikelihoodPIDAlg.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

// ROOT includes
#include "TFile.h"
#include "TMath.h"
#include "TProfile.h"

// Framework includes
#include "canvas/Persistency/Common/Ptr.h"
#include "cetlib/search_path.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/GeometryUtilities.h"

#include "cetlib/pow.h"

//------------------------------------------------------------------------------
pid::LikelihoodPIDAlg::LikelihoodPIDAlg(fhicl::ParameterSet const& pset)
{
  fmaxrr = pset.get<float>("maxrr");

  map_PhysdEdx[13] = new PhysdEdx(13);     // == muon
  map_PhysdEdx[211] = new PhysdEdx(211);   // == charged pion
  map_PhysdEdx[2212] = new PhysdEdx(2212); // == proton
}

//------------------------------------------------------------------------------
std::bitset<8> pid::LikelihoodPIDAlg::GetBitset(geo::PlaneID planeID)
{

  std::bitset<8> thisBitset;

  thisBitset.set(planeID.Plane);

  return thisBitset;
}

//------------------------------------------------------------------------------
anab::ParticleID pid::LikelihoodPIDAlg::DoParticleID(
  const std::vector<art::Ptr<anab::Calorimetry>>& calos)
{

  std::vector<anab::sParticleIDAlgScores> AlgScoresVec;
  geo::PlaneID plid;

  for (size_t i_calo = 0; i_calo < calos.size(); i_calo++) {

    art::Ptr<anab::Calorimetry> calo = calos.at(i_calo);
    if (i_calo == 0)
      plid = calo->PlaneID();
    else if (plid != calo->PlaneID())
      throw cet::exception("LikelihoodPIDAlg")
        << "PlaneID mismatch: " << plid << ", " << calo->PlaneID();

    int nptmu = 0;
    int nptpi = 0;
    int nptpro = 0;

    double lambdamu = 0;
    double lambdapi = 0;
    double lambdapro = 0;

    std::vector<float> trkdedx = calo->dEdx();
    std::vector<float> trkres = calo->ResidualRange();
    std::vector<float> trkpitches = calo->TrkPitchVec();

    for (unsigned i = 0; i < trkdedx.size(); ++i) { //hits
      //ignore the first and the last point
      if (i == 0 || i == trkdedx.size() - 1) continue;

      float this_dedx = trkdedx[i];
      float this_res = trkres[i];
      float this_pitch = trkpitches[i];

      if (this_res > fmaxrr) continue;

      // in MeV/c unit
      float p_mu = 1000. * tmc.GetTrackMomentum(this_res, 13);
      float p_pi = 1000. * tmc.GetTrackMomentum(this_res, 211);
      float p_pro = 1000. * tmc.GetTrackMomentum(this_res, 2212);

      float ke_mu = map_PhysdEdx[13]->MomentumtoKE(p_mu);
      float ke_pi = map_PhysdEdx[13]->MomentumtoKE(p_pi);
      float ke_pro = map_PhysdEdx[13]->MomentumtoKE(p_pro);

      // == muon likelihoods
      double mu_pdf, mu_pdf_max;
      bool mu_valid_pdf =
        map_PhysdEdx[13]->dEdx_PDF(ke_mu, this_pitch, this_dedx, &mu_pdf, &mu_pdf_max);
      if (mu_valid_pdf) {
        double this_lambdamu = 2. * (log(mu_pdf_max) - log(mu_pdf));
        lambdamu += this_lambdamu;
        ++nptmu;
      }

      // == pion likelihoods
      double pi_pdf, pi_pdf_max;
      bool pi_valid_pdf =
        map_PhysdEdx[211]->dEdx_PDF(ke_pi, this_pitch, this_dedx, &pi_pdf, &pi_pdf_max);
      if (pi_valid_pdf) {
        double this_lambdapi = 2. * (log(pi_pdf_max) - log(pi_pdf));
        lambdapi += this_lambdapi;
        ++nptpi;
      }

      // == proton likelihoods
      double pro_pdf, pro_pdf_max;
      bool pro_valid_pdf =
        map_PhysdEdx[2212]->dEdx_PDF(ke_pro, this_pitch, this_dedx, &pro_pdf, &pro_pdf_max);
      if (pro_valid_pdf) {
        double this_lambdapro = 2. * (log(pro_pdf_max) - log(pro_pdf));
        lambdapro += this_lambdapro;
        ++nptpro;
      }
    }

    anab::sParticleIDAlgScores lambdamuon;
    anab::sParticleIDAlgScores lambdapion;
    anab::sParticleIDAlgScores lambdaproton;

    if (nptmu) {
      lambdamuon.fAlgName = "Likelihood";
      lambdamuon.fVariableType = anab::kGOF;
      lambdamuon.fTrackDir = anab::kForward;
      lambdamuon.fAssumedPdg = 13;
      lambdamuon.fPlaneMask = GetBitset(calo->PlaneID());
      lambdamuon.fNdf = nptmu;
      lambdamuon.fValue = lambdamu / nptmu;

      AlgScoresVec.push_back(lambdamuon);
    }

    if (nptpi) {
      lambdapion.fAlgName = "Likelihood";
      lambdapion.fVariableType = anab::kGOF;
      lambdapion.fTrackDir = anab::kForward;
      lambdapion.fAssumedPdg = 221;
      lambdapion.fPlaneMask = GetBitset(calo->PlaneID());
      lambdapion.fNdf = nptpi;
      lambdapion.fValue = lambdapi / nptpi;

      AlgScoresVec.push_back(lambdapion);
    }

    if (nptpro) {
      lambdaproton.fAlgName = "Likelihood";
      lambdaproton.fVariableType = anab::kGOF;
      lambdaproton.fTrackDir = anab::kForward;
      lambdaproton.fAssumedPdg = 2212;
      lambdaproton.fPlaneMask = GetBitset(calo->PlaneID());
      lambdaproton.fNdf = nptpro;
      lambdaproton.fValue = lambdapro / nptpro;

      AlgScoresVec.push_back(lambdaproton);
    }
  }

  anab::ParticleID pidOut(AlgScoresVec, plid);

  return pidOut;
}
