/**
 * \file RiseTimeThreshold_tool.cc
 *
 * \brief Rise time is defined as the time slot in which the pulse
 * goes above a certain fraction of the maximum ADC peak value
 * given by the "PeakRatio" fhicl parameter
 *
 * @author Fran Nicolas, June 2022
 */

#include "art/Utilities/ToolConfigTable.h"
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/types/Atom.h"

#include "RiseTimeCalculatorBase.h"

namespace pmtana {

  class RiseTimeThreshold : RiseTimeCalculatorBase {

  public:
    //Configuration parameters
    struct Config {
      fhicl::Atom<double> PeakRatio{fhicl::Name("PeakRatio")};
      fhicl::Atom<bool> InterpolateSamples{fhicl::Name("InterpolateSamples"), false};
    };

    // Default constructor
    explicit RiseTimeThreshold(art::ToolConfigTable<Config> const& config);

    // Method to calculate the OpFlash t0
    double RiseTime(const pmtana::Waveform_t& wf_pulse,
                    const pmtana::PedestalMean_t& ped_pulse,
                    bool _positive) const override;

  private:
    double InterpolateTicks(size_t i, double yi, double yii, double thr) const;
    double fPeakRatio;
    bool fInterpolateSample;
  };

  RiseTimeThreshold::RiseTimeThreshold(art::ToolConfigTable<Config> const& config)
    : fPeakRatio{config().PeakRatio()},
      fInterpolateSample{config().InterpolateSamples()}
  {}

  double RiseTimeThreshold::RiseTime(const pmtana::Waveform_t& wf_pulse,
                                     const pmtana::PedestalMean_t& ped_pulse,
                                     bool _positive) const
  {

    // Pedestal-subtracted pulse
    std::vector<double> wf_aux(ped_pulse);
    if (_positive) {
      for (size_t ix = 0; ix < wf_aux.size(); ix++) {
        wf_aux[ix] = ((double)wf_pulse[ix]) - wf_aux[ix];
      }
    }
    else {
      for (size_t ix = 0; ix < wf_aux.size(); ix++) {
        wf_aux[ix] = wf_aux[ix] - ((double)wf_pulse[ix]);
      }
    }

    // rise is first sample after threshold is crossed
    auto it_max = max_element(wf_aux.begin(), wf_aux.end());
    double rise = std::lower_bound(wf_aux.begin(), it_max, fPeakRatio * (*it_max)) - wf_aux.begin();

    // linear interpolation
    if(fInterpolateSample && rise > 0){
      rise = InterpolateTicks(rise-1, wf_aux[rise-1], wf_aux[rise], fPeakRatio * (*it_max) );
    }

    return rise;
  }

  double RiseTimeThreshold::InterpolateTicks(size_t i, double y1, double y2, double thr) const 
  {
    // Linear interpolation to find x at y=thr
    // between (i,y1) and (i+1,y2)
    double frac = (thr-y1)/(y2-y1);
    return i + frac;
  }

}

DEFINE_ART_CLASS_TOOL(pmtana::RiseTimeThreshold)
