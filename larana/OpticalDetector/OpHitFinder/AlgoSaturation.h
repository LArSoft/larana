////////////////////////////////////////////////////////////////////////
// AlgoSaturation.h
// This is a hit finding algorithm adapted from AlgoSiPM.h
// This algorithm adds a correction to saturated peaks.
// The specific correction is left to be implemented by the user
// using the SaturationCorrectionTool.
//
// Created: June 16 by Laura Paulucci (lpaulucc@fnal.gov)
////////////////////////////////////////////////////////////////////////

#ifndef ALGOSATURATION_H
#define ALGOSATURATION_H

namespace fhicl {
  class ParameterSet;
}

#include "PMTPulseRecoBase.h"
#include "SaturationCorrectionTool.h"
#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"
#include "larana/OpticalDetector/OpHitFinder/SaturationCorrectionTool.h"

#include <map>
#include <memory>
#include <string>

namespace pmtana {

  class AlgoSaturation : public PMTPulseRecoBase {

  public:
    AlgoSaturation(const fhicl::ParameterSet& pset,
                   std::unique_ptr<pmtana::RiseTimeCalculatorBase> risetimecalculator = nullptr,
                   const std::string name = "AlgoSaturation");

    // Implementation of PMTPulseRecoBase::Reset() method
    void Reset();

    // Retrieving the channel associated with the pulse
    void GetWvfChannel(size_t ch) { _pulse_ch = ch; }

  protected:
    bool RecoPulse(const pmtana::Waveform_t&,
                   const pmtana::PedestalMean_t&,
                   const pmtana::PedestalSigma_t&);

    // A variable holder for a user-defined absolute ADC threshold value
    double _adc_thres;

    // Minimum width for a hit to be recorded
    int _min_width;

    // Start recording hit information after this threshold is reached
    double _2nd_thres;

    // Use this pedestal instead of the one given by the pedestal algorithm
    double _pedestal;

  private:
    std::unique_ptr<pmtana::SaturationCorrectionTool> fCorrectionTool;
    size_t _pulse_ch; //Variable to store the pulse channel
    std::vector<int>
      _ignoreChannels; //channels to which the saturation correction will not be applied
  };

}

#endif
