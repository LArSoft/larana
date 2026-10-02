// -*- mode: c++; c-basic-offset: 2; -*-
/*!
 * Title:   OpHit Algorithims
 * Authors, editors:  Ben Jones, MIT
 *                    Wes Ketchum wketchum@lanl.gov
 *                    Gleb Sinev  gleb.sinev@duke.edu
 *                    Alex Himmel ahimmel@fnal.gov
 *                    Kevin Wood  kevin.wood@stonybrook.edu
 *
 * Description:
 * These are the algorithms used by OpHitFinder to produce optical hits.
 */

#include "OpHitAlg.h"

#include "larana/OpticalDetector/OpHitFinder/PulseRecoManager.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larreco/Calibrator/IPhotonCalibrator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <vector>

namespace opdet {
  //----------------------------------------------------------------------------
  void RunHitFinder(std::vector<raw::OpDetWaveform> const& opDetWaveformVector,
                    std::vector<recob::OpHit>& hitVector,
                    pmtana::PulseRecoManager const& pulseRecoMgr,
                    pmtana::PMTPulseRecoBase const& threshAlg,
                    geo::WireReadoutGeom const& wireReadoutGeom,
                    float hitThreshold,
                    detinfo::DetectorClocksData const& clocksData,
                    calib::IPhotonCalibrator const& calibrator,
                    bool use_start_time,
                    bool timestamp_is_relative)
  {
    for (auto const& waveform : opDetWaveformVector) {
      const int channel = static_cast<int>(waveform.ChannelNumber());

      if (!wireReadoutGeom.IsValidOpChannel(channel)) {
        mf::LogError("OpHitFinder")
          << "Error! unrecognized channel number " << channel << ". Ignoring pulse";
        continue;
      }

      pulseRecoMgr.Reconstruct(waveform);

      // Get the result
      auto const& pulses = threshAlg.GetPulses();

      const double timeStamp = waveform.TimeStamp();

      for (auto const& pulse : pulses)
        ConstructHit(hitThreshold,
                     channel,
                     timeStamp,
                     pulse,
                     hitVector,
                     clocksData,
                     calibrator,
                     use_start_time,
                     timestamp_is_relative);
    }
  }

  //----------------------------------------------------------------------------
  void ConstructHit(float hitThreshold,
                    int channel,
                    double timeStamp,
                    pmtana::pulse_param const& pulse,
                    std::vector<recob::OpHit>& hitVector,
                    detinfo::DetectorClocksData const& clocksData,
                    calib::IPhotonCalibrator const& calibrator,
                    bool use_start_time,
                    bool timestamp_is_relative)
  {
    if (pulse.peak < hitThreshold) return;

    auto tick_period = clocksData.OpticalClock().TickPeriod();
    double absTime = timeStamp +  tick_period * (use_start_time ? pulse.t_start : pulse.t_max);
    double startTime = timeStamp + tick_period * pulse.t_start;
    // ProtoDUNE-HD/VD OpDetWaveform timestamps were made relative to the trigger timestamp in the following commit
    // https://github.com/DUNE/duneprototypes/pull/109/changes/6ab18a41cdc78edcc609bc5043ad93d9b5ac2134
    // So we can turn off this shifting here if necessary
    if (timestamp_is_relative) {
       absTime += clocksData.TriggerTime();
    }
    else {
       startTime -= clocksData.TriggerTime();
    }
    double relTime = absTime - clocksData.TriggerTime();

    double riseTime = tick_period * pulse.t_rise;

    int frame = clocksData.OpticalClock().Frame(timeStamp);

    double PE = 0.0;
    if (calibrator.UseArea())
      PE = calibrator.PE(pulse.area, channel);
    else
      PE = calibrator.PE(pulse.peak, channel);

    double width = (pulse.t_end - pulse.t_start) * tick_period;

    hitVector.emplace_back(channel,
                           relTime,
                           absTime,
                           startTime,
                           riseTime,
                           frame,
                           width,
                           pulse.area,
                           pulse.peak,
                           PE,
                           0.0);
  }

} // End namespace opdet
