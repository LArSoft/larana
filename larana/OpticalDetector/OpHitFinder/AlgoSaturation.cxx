/////////////////////////////////////////////////////////////////////////
// AlgoSaturation.cxx
//
// Created: Aug 17 by Laura Paulucci (lpaulucc@fnal.gov)
//
// This algorithm identifies OpticalHits and corrects saturation using a ToS method.
////////////////////////////////////////////////////////////////////////

#include "AlgoSaturation.h"

#include "fhiclcpp/ParameterSet.h"

namespace pmtana {

  //---------------------------------------------------------------------------
  AlgoSaturation::AlgoSaturation(const fhicl::ParameterSet& pset,
                                 std::unique_ptr<pmtana::RiseTimeCalculatorBase> risetimecalculator,
                                 const std::string name)
    : PMTPulseRecoBase(name)
  {

    _adc_thres = pset.get<float>("ADCThreshold");
    _min_width = pset.get<float>("MinWidth");
    _2nd_thres = pset.get<float>("SecondThreshold");
    _pedestal = pset.get<float>("Pedestal");
    _ignoreChannels = pset.get<std::vector<int>>("IgnoreChannels");
    _calibThrs = pset.get<double>("calibThrs"); //Value at which the calibration curve was evaluated
    _channelVec = pset.get<std::vector<size_t>>("channelVec");
    _calibParam1 = pset.get<std::vector<double>>("calibParam1");
    _calibParam2 = pset.get<std::vector<double>>("calibParam2");
    _calibParam3 = pset.get<std::vector<double>>("calibParam3");

    _risetime_calc_ptr = std::move(risetimecalculator);
    _min_plateau_size = pset.get<int>("MinPlateauSize", 2); // Min consecutive flat samples

    Reset();
  }

  //---------------------------------------------------------------------------
  void AlgoSaturation::Reset()
  {
    PMTPulseRecoBase::Reset();
  }

  //---------------------------------------------------------------------------
  bool AlgoSaturation::RecoPulse(const pmtana::Waveform_t& wf,
                                 const pmtana::PedestalMean_t& ped_mean,
                                 const pmtana::PedestalSigma_t& ped_rms)
  {

    bool fire = false;
    bool record_hit = false;
    int counter = 0;
    double pedestal =
      ped_mean
        .front(); //Switch pedestal definition to incoroprate pedestal finder - K.S. 04/18/2019

    double threshold = _adc_thres;
    threshold += pedestal;
    double pre_threshold = _2nd_thres;
    pre_threshold += pedestal;

    //Reset();
    _num_plateaus.clear();
    _plateau_lengths.clear();

    //Create the calibration map
    for (size_t i = 0; i < _channelVec.size(); i++) {
      fcalibMap[_channelVec.at(i)] = {_calibParam1.at(i), _calibParam2.at(i), _calibParam3.at(i)};
    }

    for (short const& value : wf) {

      // Retrieve the channel number set prior to running reconstruction
      const int ch = Channel();

      if (!fire && (double(value) >= pre_threshold)) {

        // Found a new pulse
        fire = true;
        record_hit = false;
        _pulse.t_start = counter;
      }

      if (fire && (double(value) < pre_threshold)) {

        // Found the end of a pulse
        fire = false;
        _pulse.t_end = counter - 1;
        if (record_hit && ((_pulse.t_end - _pulse.t_start) >= _min_width)) {
          //Check for saturation
          if (std::find(_ignoreChannels.begin(), _ignoreChannels.end(), ch) ==
              _ignoreChannels.end()) {
            CorrectSaturation(wf, _pulse, ch);
          }
          if (_risetime_calc_ptr)
            _pulse.t_rise = _risetime_calc_ptr->RiseTime(
              {wf.begin() + _pulse.t_start, wf.begin() + _pulse.t_end},
              {ped_mean.begin() + _pulse.t_start, ped_mean.begin() + _pulse.t_end},
              true);

          _pulse_v.push_back(_pulse); //PROBLEMAS AQUI!!!
          record_hit = false;
        }
        _pulse.reset_param();
      }

      if (fire) {
        // We want to record the hit only if _adc_thres is reached
        if (!record_hit && (double(value) >= threshold)) record_hit = true;

        // Add this ADC count to the integral
        _pulse.area += (double(value) - double(pedestal));

        if (_pulse.peak < (double(value) - double(pedestal))) {

          // Found a new maximum
          _pulse.peak = (double(value) - double(pedestal));
          _pulse.t_max = counter;
        }
      }
      counter++;
    }

    if (fire) {

      // Take care of a pulse that did not finish within the readout window
      fire = false;
      _pulse.t_end = counter - 1;
      if (record_hit && ((_pulse.t_end - _pulse.t_start) >= _min_width)) {
        if (_risetime_calc_ptr)
          _pulse.t_rise = _risetime_calc_ptr->RiseTime(
            {wf.begin() + _pulse.t_start, wf.begin() + _pulse.t_end},
            {ped_mean.begin() + _pulse.t_start, ped_mean.begin() + _pulse.t_end},
            true);

        _pulse_v.push_back(_pulse);
        record_hit = false;
      }
      _pulse.reset_param();
    }

    return true;
  }

  void AlgoSaturation::CorrectSaturation(pmtana::Waveform_t const& wf,
                                         pmtana::pulse_param& _pulse,
                                         int ch)
  {
    bool in_saturation = false;
    int satcounter = 0, plateau_count = 0;
    std::vector<int> plateau_lengths;
    double tolerance = 0.1; //to be optimized

    // Ensure valid indices
    if (_pulse.t_start >= _pulse.t_end || _pulse.t_end >= static_cast<int>(wf.size())) return;
    // Start loop at max(1, t_start) to prevent wf[i-1] out-of-bounds access
    size_t start_tick = std::max(static_cast<size_t>(1), static_cast<size_t>(_pulse.t_start));

    for (size_t i = start_tick; i <= static_cast<size_t>(_pulse.t_end);
         ++i) { //go over the full lenght of a pulse looking for saturation (slope ~ 0 near peak)
      if ((static_cast<double>(wf[i]) - static_cast<double>(wf[i - 1]) <= tolerance) &&
          (static_cast<double>(wf[i]) >=
           _calibThrs)) { //the last one is just a safegard against spurious fluctuations
        if (!in_saturation) {
          in_saturation = true;
          satcounter = 1;
        }
        else {
          satcounter++;
        }
      }
      else { // End of a saturated region
        if (in_saturation) {
          if (satcounter >= _min_plateau_size) {
            plateau_count++;
            plateau_lengths.push_back(satcounter);
          }
          in_saturation = false;
          satcounter = 0;
        }
      }
    }

    // Catch trailing saturation at the end of pulse window
    if (in_saturation && satcounter >= _min_plateau_size) { plateau_lengths.push_back(satcounter); }

    if (!plateau_lengths.empty()) {
      //For the moment, the correction applied here is based on the "time-over-saturation (ToS)" method, where the missing area is added to the original area calculation
      //(see Maressa's presentation here: https://indico.fnal.gov/event/74741/contributions/342609/attachments/198490/276375/Jun_ToTUpdate_PDVD.pdf)
      //parameters of the calib curve, channel dependent, in order: a, b, c
      std::vector<double> calibvec = fcalibMap[ch];
      for (size_t i = 0; i < plateau_lengths.size(); i++) {
        double missingArea = calibvec[0] * pow(plateau_lengths[i], calibvec[1]) *
                             exp(calibvec[2] * plateau_lengths[i]);
        //for debugging
        std::cout << "Saturated peak " << i << " ch " << ch << " old area " << _pulse.area
                  << " sat ticks " << plateau_lengths[i] << " missing area " << missingArea
                  << std::endl;
        _pulse.area += missingArea;
      }
    }
  }

}
