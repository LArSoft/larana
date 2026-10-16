/////////////////////////////////////////////////////////////////////////
// AlgoSaturation.cxx
//
// Created: Aug 17 by Laura Paulucci (lpaulucc@fnal.gov)
//
// This algorithm identifies OpticalHits and allow for the correction
// of saturated peaks.
////////////////////////////////////////////////////////////////////////

#include "AlgoSaturation.h"
#include "art/Utilities/make_tool.h"
#include "fhiclcpp/ParameterSet.h"
#include <algorithm>

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

    _risetime_calc_ptr = std::move(risetimecalculator);

    // Dynamically load whichever art::Tool plugin is specified in FHiCL
    if (pset.has_key("CorrectionTool")) {
      const auto tool_pset = pset.get<fhicl::ParameterSet>("CorrectionTool");
      fCorrectionTool = art::make_tool<pmtana::SaturationCorrectionTool>(tool_pset);
    }
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
          //Check for saturation: Execute tool if configured and channel is not ignored
          if (fCorrectionTool && std::find(_ignoreChannels.begin(), _ignoreChannels.end(), ch) ==
                                   _ignoreChannels.end()) {
            fCorrectionTool->CorrectSaturation(wf, _pulse, ch);
          }
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

}
