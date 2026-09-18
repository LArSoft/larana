/**
 * @file   larana/OpticalDetector/AlgoSaturationMaker_tool.cc
 * @brief  _art_ tool to create a `pmtana::AlgoSaturation` algorithm.
 * @author Laura Paulucci (lpaulucc@fnal.gov)
 * @date   September 3, 2026
 */

// LArSoft libraries
#include "larana/OpticalDetector/HitAlgoMakerToolBase.h"
#include "larana/OpticalDetector/OpHitFinder/AlgoSaturation.h"

// framework libraries
#include "art/Utilities/ToolMacros.h"

// -----------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(opdet::HitAlgoMakerToolBase<pmtana::AlgoSaturation>)
