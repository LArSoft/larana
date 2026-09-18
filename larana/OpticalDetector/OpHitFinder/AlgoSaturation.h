////////////////////////////////////////////////////////////////////////
// AlgoSaturation.h
// This is a hit finding algorithm adapted from AlgoSiPM.h
// The original code makes an optical hit out of
// everything above threshold, but uses only the first peak to assign hit time.
// This algorithm adds a correction to saturated peaks using a ToS method.
//
// Created: June 16 by Laura Paulucci (lpaulucc@fnal.gov)
////////////////////////////////////////////////////////////////////////

#ifndef ALGOSATURATION_H
#define ALGOSATURATION_H

namespace fhicl {
  class ParameterSet;
}

#include "PMTPulseRecoBase.h"
#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"

#include <string>
#include <map>

namespace pmtana {

  class AlgoSaturation : public PMTPulseRecoBase {

  public:
    AlgoSaturation(const fhicl::ParameterSet& pset,
             std::unique_ptr<pmtana::RiseTimeCalculatorBase> risetimecalculator = nullptr,
             const std::string name = "AlgoSaturation");

    // Implementation of PMTPulseRecoBase::Reset() method
    void Reset();

    // Retrieving the channel associated with the pulse
    void GetWvfChannel(size_t ch) {_pulse_ch = ch;}

    // A method to set user-defined ADC threshold value
    //      void SetADCThreshold(double v) {_adc_thres = v;};

    // A method to set a multiplication factor to the pedestal standard deviation
    // which is used as one of two input values to define a threshold
    //      void SetNSigma(double v) {_nsigma = v;};

    // Methods for saturation results
    const std::vector<std::vector<int>>& GetPlateauLengths() const { return _plateau_lengths; }
    const std::vector<int>& GetNumPlateaus() const { return _num_plateaus; }

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

    // A variable holder for a multiplicative factor for the pedestal
    // standard deviation to define the threshold
    //      double _nsigma;
    void CorrectSaturation(pmtana::Waveform_t const& wf, pmtana::pulse_param& _pulse, int ch);

  private:
    size_t _pulse_ch; //Variable to store the pulse channel  
 //   double _sat_threshold; // Threshold value (ADC or Amplitude) defining saturation
    int _min_plateau_size; // Minimum consecutive ADC samples to count as a plateau
    // Stores plateau statistics per reconstructed hit
    std::vector<int> _num_plateaus;                // Number of plateaus in each hit
    std::vector<std::vector<int>> _plateau_lengths; // Duration (samples) of each plateau per hit
    std::vector<int> _ignoreChannels; //channels to which the saturation correction will not be applied
    double _calibThrs;     //Value at which the calibration curve was evaluated
    //list of channels and calibration parameters
    std::vector<size_t> _channelVec;
    std::vector<double> _calibParam1;
    std::vector<double> _calibParam2;
    std::vector<double> _calibParam3;
    std::map<size_t,std::vector<double>> fcalibMap;
  };

}

#endif
