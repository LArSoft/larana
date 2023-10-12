////////////////////////////////////////////////////////////////////////
// \file OpDetResponse.h
//
// \brief service containing information about the response of optical detectors
//
// \author ahimmel@phy.duke.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef OPDET_RESPONSE_INTERFACE_H
#define OPDET_RESPONSE_INTERFACE_H

// LArSoft includes
#include "larcore/Geometry/WireReadout.h"
namespace sim {
  struct OnePhoton;
}

// ART includes
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/fwd.h"

namespace opdet {
  class OpDetResponseInterface {
  public:
    virtual ~OpDetResponseInterface() = default;

    virtual void reconfigure(fhicl::ParameterSet const& p);

    virtual int NOpChannels() const;

    // Can only uniquely go from readout to geometry since
    // one geometrical channel goes to multiple readout channels
    virtual int readoutToGeoChannel(int readoutChannel) const;

    virtual bool detected(int OpChannel, const sim::OnePhoton& Phot, int& newOpChannel) const;
    virtual bool detected(int OpChannel, const sim::OnePhoton& Phot) const;
    virtual bool detectedLite(int OpChannel, int& newOpChannel) const;
    virtual bool detectedLite(int OpChannel) const;

    virtual float wavelength(double energy) const;

  private:
    virtual void doReconfigure(fhicl::ParameterSet const& p) = 0;

    virtual int doNOpChannels() const;
    virtual int doReadoutToGeoChannel(int readoutChannel) const;

    virtual bool doDetected(int OpChannel, const sim::OnePhoton& Phot, int& newOpChannel) const = 0;
    virtual bool doDetectedLite(int OpChannel, int& newOpChannel) const = 0;

  }; // class OpDetResponse

  //-------------------------------------------------------------------------------------------------------------
  inline void OpDetResponseInterface::reconfigure(fhicl::ParameterSet const& p)
  {
    doReconfigure(p);
  }

  //-------------------------------------------------------------------------------------------------------------
  inline int OpDetResponseInterface::NOpChannels() const { return doNOpChannels(); }

  //-------------------------------------------------------------------------------------------------------------
  inline int OpDetResponseInterface::doNOpChannels() const
  {
    // By default return the number of detector channels
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();
    return wireReadoutGeom.NOpChannels();
  }

  //-------------------------------------------------------------------------------------------------------------
  inline int OpDetResponseInterface::readoutToGeoChannel(int readoutChannel) const
  {
    return doReadoutToGeoChannel(readoutChannel);
  }

  //-------------------------------------------------------------------------------------------------------------
  inline int OpDetResponseInterface::doReadoutToGeoChannel(int readoutChannel) const
  {
    auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();
    return wireReadoutGeom.OpDetFromOpChannel(readoutChannel);
  }

  //-------------------------------------------------------------------------------------------------------------
  inline bool OpDetResponseInterface::detected(int OpChannel,
                                               const sim::OnePhoton& Phot,
                                               int& newOpChannel) const
  {
    return doDetected(OpChannel, Phot, newOpChannel);
  }

  //-------------------------------------------------------------------------------------------------------------
  inline bool OpDetResponseInterface::detected(int OpChannel, const sim::OnePhoton& Phot) const
  {
    int newOpChannel;
    return doDetected(OpChannel, Phot, newOpChannel);
  }

  //-------------------------------------------------------------------------------------------------------------
  inline bool OpDetResponseInterface::detectedLite(int OpChannel, int& newOpChannel) const
  {
    return doDetectedLite(OpChannel, newOpChannel);
  }

  //-------------------------------------------------------------------------------------------------------------
  inline bool OpDetResponseInterface::detectedLite(int OpChannel) const
  {
    int newOpChannel;
    return doDetectedLite(OpChannel, newOpChannel);
  }

  //-------------------------------------------------------------------------------------------------------------
  inline float OpDetResponseInterface::wavelength(double energy) const
  {
    return (2.0 * 3.142) * 0.000197 / energy;
  }

} //namespace opdet

DECLARE_ART_SERVICE_INTERFACE(opdet::OpDetResponseInterface, LEGACY)

#endif //OPDET_RESPONSE_H
