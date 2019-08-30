cmake_minimum_required(VERSION 2.8.7)

get_filename_component(dirName ${CMAKE_CURRENT_LIST_DIR} NAME)
if(thisTarget)
  set (thisTarget "${thisTarget}.${dirName}")
else(thisTarget)
  set (thisTarget "${dirName}")
endif(thisTarget)
message(STATUS "Source : Specified  /  Target : ${thisTarget}")

add_library(${thisTarget} OBJECT
  NeuKine/ccapnu.f
  NeuKine/cmucap.f
  NeuKine/cpimuPolari.f
  NeuKine/csNeumuEMu.f
  NeuKine/csampNeueEMu.f
  cDeltaDecay.f
  cDsDecay.f
  cbomegaDcy.f
  cdDecay.f
  cetaDecay.f
  cetapDecay.f
  cgzaiDecay.f
  ckaonDecay.f
  clambdaDcy.f
  clambdacDcy.f
  cmuNeuDcy.f
  cpCos2pxyz.f
  cpi0Decay.f
  cpiMuDecay.f
  csigmaDecay.f
  ctauNeuDcy.f
  cvmDecay.f
  )
