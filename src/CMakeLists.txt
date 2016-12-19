cmake_minimum_required( VERSION 3.1.0 FATAL_ERROR )

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

find_package(ROOT REQUIRED Core Matrix Physics MathCore SPlot MathMore Graf RooFit RooStats RIO Gpad Hist Tree TreePlayer TMVA XrdUtils NetxNG Postscript ASImage HistPainter)
include_directories( ${ROOT_INCLUDE_DIR} )
if (ROOT_FOUND)
    MESSAGE("--   ROOT include directories are " ${ROOT_INCLUDE_DIRS})
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -stdlib=libc++")
endif()

ROOT_GENERATE_DICTIONARY( RooMisID_dict RooMomentumFractionPdf.h RooMisIDGaussian.h RooPowerLawApproxVar.h RooInverseGaussian.h )
add_library( RooMisID SHARED RooMomentumFractionPdf.cxx RooMisIDGaussian.cxx RooPowerLawApproxVar.cxx RooInverseGaussian.cxx RooMisID_dict.cxx )
target_link_libraries( RooMisID ${ROOT_LIBRARIES} RooFitCore )
