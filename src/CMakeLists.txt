cmake_minimum_required( VERSION 2.8 )

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

find_package( ROOT REQUIRED )
include_directories( ${ROOT_INCLUDE_DIR} )

ROOT_GENERATE_DICTIONARY( RooMisID_dict RooMomentumFractionPdf.h RooMisIDGaussian.h RooPowerLawApproxVar.h RooInverseGaussian.h )
add_library( RooMisID SHARED RooMomentumFractionPdf.cxx RooMisIDGaussian.cxx RooPowerLawApproxVar.cxx RooInverseGaussian.cxx RooMisID_dict.cxx )
target_link_libraries( RooMisID ${ROOT_LIBRARIES} RooFitCore )
