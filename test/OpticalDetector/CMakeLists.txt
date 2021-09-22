# ======================================================================
#
# Testing
#
# ======================================================================

include(CetTest)
cet_enable_asserts()

# test needs fixing after the OpticalRecoUpdate changes
#cet_test(AlgoThreshold_test USE_BOOST_UNIT
#			    LIBRARIES larana_OpticalDetector
#					     fhiclcpp::fhiclcpp
#)

cet_test(OpFlashAlg_test USE_BOOST_UNIT
			 LIBRARIES larana_OpticalDetector
)

#cet_test(standalone_test)
