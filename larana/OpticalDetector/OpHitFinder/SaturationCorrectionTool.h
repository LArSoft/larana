/////////////////////////////////////////////////////////////////////////
// SaturationCorrectionTool.h
//
// Created: Oct 14 by Laura Paulucci (lpaulucc@fnal.gov)
//
// This is a base tool for saturation correction.
// The user will select which algorithm to run (to be
// implemented in the specific experiment repository)
// directly from the FHiCL configuration file.
// This allows for different saturation correction methods to be implemented.
////////////////////////////////////////////////////////////////////////

#ifndef SATURATIONCORRECTIONTOOL_H
#define SATURATIONCORRECTIONTOOL_H

#include "larana/OpticalDetector/OpHitFinder/OpticalRecoTypes.h"
#include "larana/OpticalDetector/OpHitFinder/PMTPulseRecoBase.h"

namespace pmtana {

  class SaturationCorrectionTool {
  public:
    virtual ~SaturationCorrectionTool() = default;

    // Pure virtual interface method to override in experiment-specific tools
    virtual void CorrectSaturation(pmtana::Waveform_t const& wf,
                                   pmtana::pulse_param& pulse,
                                   double pedestal,
                                   int ch) = 0;
  };

} // namespace pmtana

#endif
