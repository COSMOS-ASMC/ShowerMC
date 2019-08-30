### CMake project file for COSMOS
cmake_minimum_required (VERSION 2.8.7)


### set default compiler as intel compiler for MPI if exist
if(NOT CMAKE_Fortran_COMPILER AND NOT CMAKE_Fortran_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_Fortran_COMPILER NAMES mpiifort DOC "Fortran compiler")
  if(CMAKE_Fortran_COMPILER)
    message(STATUS "Intel Fortran Compiler ${CMAKE_Fortran_COMPILER} is found")
  endif(CMAKE_Fortran_COMPILER)
endif(NOT CMAKE_Fortran_COMPILER AND NOT CMAKE_Fortran_COMPILER_INIT)
if(NOT CMAKE_C_COMPILER AND NOT CMAKE_C_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_C_COMPILER NAMES mpiicc DOC "C compiler")
  if(CMAKE_C_COMPILER)
    message(STATUS "Intel C Compiler ${CMAKE_C_COMPILER} is found")
  endif(CMAKE_C_COMPILER)
endif(NOT CMAKE_C_COMPILER AND NOT CMAKE_C_COMPILER_INIT)
if(NOT CMAKE_CXX_COMPILER AND NOT CMAKE_CXX_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_CXX_COMPILER NAMES mpiicpc DOC "CXX compiler")
  if(CMAKE_CXX_COMPILER)
    message(STATUS "Intel CXX Compiler ${CMAKE_CXX_COMPILER} is found")
  endif(CMAKE_CXX_COMPILER)
endif(NOT CMAKE_CXX_COMPILER AND NOT CMAKE_CXX_COMPILER_INIT)

### set default compiler as intel compiler if exist
if(NOT CMAKE_Fortran_COMPILER AND NOT CMAKE_Fortran_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_Fortran_COMPILER NAMES ifort DOC "Fortran compiler")
  if(CMAKE_Fortran_COMPILER)
    message(STATUS "Intel Fortran Compiler ${CMAKE_Fortran_COMPILER} is found")
  endif(CMAKE_Fortran_COMPILER)
endif(NOT CMAKE_Fortran_COMPILER AND NOT CMAKE_Fortran_COMPILER_INIT)
if(NOT CMAKE_C_COMPILER AND NOT CMAKE_C_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_C_COMPILER NAMES icc DOC "C compiler")
  if(CMAKE_C_COMPILER)
    message(STATUS "Intel C Compiler ${CMAKE_C_COMPILER} is found")
  endif(CMAKE_C_COMPILER)
endif(NOT CMAKE_C_COMPILER AND NOT CMAKE_C_COMPILER_INIT)
if(NOT CMAKE_CXX_COMPILER AND NOT CMAKE_CXX_COMPILER_INIT)
  FIND_PROGRAM(CMAKE_CXX_COMPILER NAMES icpc DOC "CXX compiler")
  if(CMAKE_CXX_COMPILER)
    message(STATUS "Intel CXX Compiler ${CMAKE_CXX_COMPILER} is found")
  endif(CMAKE_CXX_COMPILER)
endif(NOT CMAKE_CXX_COMPILER AND NOT CMAKE_CXX_COMPILER_INIT)


### set project name
project (SHOWERLIB)
enable_language (Fortran)
enable_language (C)

### set SHOWERMCTOP
set (SHOWERMCTOP $ENV{SHOWERMCTOP})
if (NOT SHOWERMCTOP)
  message(FATAL_ERROR "$SHOWERMCTOP is not defined")
endif (NOT SHOWERMCTOP)

### set COSMOSTOP AND COSMOSINC
set (COSMOSTOP $ENV{COSMOSTOP})
if (NOT COSMOSTOP)
  message(FATAL_ERROR "$COSMOSTOP is not defined")
endif (NOT COSMOSTOP)
set (COSMOSINC ${COSMOSTOP}/cosmos)

### set LIBLOFT AND LIBLOFTINC
set (LIBLOFT $ENV{LIBLOFT})
if (NOT LIBLOFT)
  message(FATAL_ERROR "$LIBLOFT is not defined")
endif (NOT LIBLOFT)
set (LIBLOFTINC ${LIBLOFT}/Header)


# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: None Debug Release."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)


# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
  "Choose the installation directory; by default it installs in the COSMOS directory."
  FORCE)


# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)


### set compile option
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (ARCH LinuxGfort)
  set (CMAKE_Fortran_FLAGS_RELEASE "-cpp -std=legacy -fno-align-commons -O3 -D${ARCH}")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-cpp -std=legacy -fno-align-commons -O0 -D${ARCH} -g")
  set (CMAKE_C_FLAGS_RELEASE "-O3 -D${ARCH}")
  set (CMAKE_C_FLAGS_DEBUG   "-O0 -D${ARCH} -g")
elseif (Fortran_COMPILER_NAME MATCHES "mpiifort.*")
  # mpiifort
  set (ARCH MPIIFC64)
  set (CMAKE_Fortran_FLAGS_RELEASE "-cpp -traceback -warn noalign -diag-disable=remark -align dcommons -O3 -D${ARCH} -DDOMPI")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-cpp -traceback -warn noalign -diag-disable=remark -align dcommons -O0 -D${ARCH} -g ")
  set (CMAKE_C_FLAGS_RELEASE "-O3 -D${ARCH}")
  set (CMAKE_C_FLAGS_DEBUG   "-O0 -D${ARCH} -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort
  find_program(MPIEXEC mpiexec_mpt)
  if (MPIEXEC)
    set (ARCH MPTIFC64)
    set (CMAKE_Fortran_FLAGS_RELEASE "-cpp -traceback -warn noalign -diag-disable=remark -align dcommons -O3 -D${ARCH} -DDOMPI")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-cpp -traceback -warn noalign -diag-disable=remark -align dcommons -O0 -D${ARCH} -g ")
    set (CMAKE_C_FLAGS_RELEASE "-O3 -D${ARCH}")
    set (CMAKE_C_FLAGS_DEBUG   "-O0 -D${ARCH} -g")
  else (MPIEXEC)
    set (ARCH PCLinuxIFC64)
    set (CMAKE_Fortran_FLAGS_RELEASE "-cpp -traceback -warn noalign -diag-disable=remark -align dcommons -O3 -D${ARCH}")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-cpp -traceback -warn noalign -diag-disable=remark -align dcommons -O0 -D${ARCH} -g ")
    set (CMAKE_C_FLAGS_RELEASE "-O3 -D${ARCH}")
    set (CMAKE_C_FLAGS_DEBUG   "-O0 -D${ARCH} -g")
  endif (MPIEXEC)
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")


### dump variables
message (STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER})
message (STATUS "C compiler: " ${CMAKE_C_COMPILER})
message (STATUS "C++ compiler: " ${CMAKE_CXX_COMPILER})
message (STATUS "$SHOWERMCTOP is ${SHOWERMCTOP}")
message (STATUS "$LIBLOFT is ${LIBLOFT}")
message (STATUS "$COSMOSTOP is ${COSMOSTOP}")
message (STATUS "$ARCH is ${ARCH}")


include_directories(${LIBLOFTINC})
include_directories(${COSMOSINC})
set (CMAKE_Fortran_MODULE_DIRECTORY "${SHOWERMCTOP}/lib/${ARCH}")
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${SHOWERMCTOP}/lib/${ARCH}")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${SHOWERMCTOP}/lib/${ARCH}")
set (LIBSHOWERMC showermc)

### function
FUNCTION (SET_INCLUDE_DIR_LIST _InFileList)
  foreach (DIR ${_InFileList})
    if(IS_DIRECTORY ${DIR})
      if ( ${DIR} MATCHES "/CMakeFiles$" OR
	   ${DIR} MATCHES "/Test$" OR
	   ${DIR} MATCHES "/Util$" OR
	   ${DIR} MATCHES "/UC$" OR
	   ${DIR} MATCHES "/Ana$" OR
	   ${DIR} MATCHES "/Ceren$" OR
	   ${DIR} MATCHES "/Data$" OR
	   ${DIR} MATCHES "/Hidden$" OR
	   ${DIR} MATCHES "/Misc$" OR
	   ${DIR} MATCHES "/MuUg$" OR
	   ${DIR} MATCHES "/Obso$" OR
	   ${DIR} MATCHES "/TeX$" OR
	   ${DIR} MATCHES "/lib$" OR
	   ${DIR} MATCHES "/UserHook$" OR
	   ${DIR} MATCHES "/Doc$" OR
	   ${DIR} MATCHES "/Betterone$" OR
	   ${DIR} MATCHES "/Scrpt$" OR
	   ${DIR} MATCHES "/NotForExport$" OR
	   ${DIR} MATCHES "/LinuxGfort$" OR
	   ${DIR} MATCHES "/${ARCH}$" 
	  )
	list (REMOVE_ITEM _InFileList ${DIR})
      else()
#	message(${DIR})
      endif()
    else(IS_DIRECTORY ${DIR})
      list (REMOVE_ITEM _InFileList ${DIR})
    endif(IS_DIRECTORY ${DIR})
  endforeach (DIR ${_InFileList})
  set(DIR_LIST ${_InFileList} PARENT_SCOPE)
ENDFUNCTION (SET_INCLUDE_DIR_LIST _InFileList)

set (SEE_SUBDIR TRUE)

get_filename_component(TARLIBLOFT ${LIBLOFT} NAME)
get_filename_component(TARCOSMOS ${COSMOSTOP} NAME)

### add subdirectory
#foreach(subdir ${LIBLOFT} )
#  add_subdirectory(${subdir})
#endforeach(subdir)

add_library(${LIBSHOWERMC}

  ### LibLoft/Decay
  ${LIBLOFT}/Decay/NeuKine/ccapnu.f
  ${LIBLOFT}/Decay/NeuKine/cmucap.f
  ${LIBLOFT}/Decay/NeuKine/cpimuPolari.f
  ${LIBLOFT}/Decay/NeuKine/csNeumuEMu.f
  ${LIBLOFT}/Decay/NeuKine/csampNeueEMu.f
  ${LIBLOFT}/Decay/cDeltaDecay.f
  ${LIBLOFT}/Decay/cDsDecay.f
  ${LIBLOFT}/Decay/cbomegaDcy.f
  ${LIBLOFT}/Decay/cdDecay.f
  ${LIBLOFT}/Decay/cetaDecay.f
  ${LIBLOFT}/Decay/cetapDecay.f
  ${LIBLOFT}/Decay/cgzaiDecay.f
  ${LIBLOFT}/Decay/ckaonDecay.f
  ${LIBLOFT}/Decay/clambdaDcy.f
  ${LIBLOFT}/Decay/clambdacDcy.f
  ${LIBLOFT}/Decay/cmuNeuDcy.f
  ${LIBLOFT}/Decay/cpCos2pxyz.f
  ${LIBLOFT}/Decay/cpi0Decay.f
  ${LIBLOFT}/Decay/cpiMuDecay.f
  ${LIBLOFT}/Decay/csigmaDecay.f
  ${LIBLOFT}/Decay/ctauNeuDcy.f
  ${LIBLOFT}/Decay/cvmDecay.f

  ### LibLoft/EM/FromEpi
  ${LIBLOFT}/EM/FromCos/MixedMCS/DCS/cDCSTypeDef.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/DCS/cDCS_Z2file.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/DCS/cDCSconstForMedia.f
  ${LIBLOFT}/EM/FromCos/MixedMCS/DCS/cManipDCS.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/DCS/cmkDCSconst.f
  ${LIBLOFT}/EM/FromCos/MixedMCS/Sampling/cMCS0.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/Sampling/cMCSTypeDef.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/Sampling/cMCSconstForMedia.f
  ${LIBLOFT}/EM/FromCos/MixedMCS/Sampling/cfixMCSmodel.f
  ${LIBLOFT}/EM/FromCos/MixedMCS/Sampling/cgetLambda.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/Sampling/cgetMuc.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/Sampling/ciniMixedMS.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/Sampling/cmkMCSconst.f
  ${LIBLOFT}/EM/FromCos/MixedMCS/Sampling/csampCS.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/Sampling/csampSoftMCS.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/TPXS/cMWconstForMedia.f
  ${LIBLOFT}/EM/FromCos/MixedMCS/TPXS/cMWscatPrep.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/TPXS/cManipTPXS.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/TPXS/cTPXSTypeDef.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/TPXS/cTPXS_Z2file.f90
  ${LIBLOFT}/EM/FromCos/MixedMCS/TPXS/cmkTPXSconst.f
  ${LIBLOFT}/EM/FromCos/Sync/cgUpsilon.f
  ${LIBLOFT}/EM/FromCos/Sync/cgetBsin.f
  ${LIBLOFT}/EM/FromCos/Sync/ck23.f
  ${LIBLOFT}/EM/FromCos/Sync/cmBremE.f
  ${LIBLOFT}/EM/FromCos/Sync/cmBremE1.f
  ${LIBLOFT}/EM/FromCos/Sync/cmBremE2.f
  ${LIBLOFT}/EM/FromCos/Sync/cmBremF1.f
  ${LIBLOFT}/EM/FromCos/Sync/cmBremF2.f
  ${LIBLOFT}/EM/FromCos/Sync/cmBremI1.f
  ${LIBLOFT}/EM/FromCos/Sync/cmBremI2.f
  ${LIBLOFT}/EM/FromCos/Sync/cmBremMFP.f
  ${LIBLOFT}/EM/FromCos/Sync/cmPairE.f
  ${LIBLOFT}/EM/FromCos/Sync/cmPairMFP.f
  ${LIBLOFT}/EM/FromCos/Sync/cmPairSpec.f
  ${LIBLOFT}/EM/FromCos/Sync/csyncTELoss.f
  ${LIBLOFT}/EM/FromCos/Sync/ctXai.f
  ${LIBLOFT}/EM/FromCos/Sync/cupsilon.f
  ${LIBLOFT}/EM/FromCos/Sync/cxai.f
  ${LIBLOFT}/EM/FromCos/cGetXXsec.f
  ${LIBLOFT}/EM/FromCos/cHowMCS.f90
  ${LIBLOFT}/EM/FromCos/cbremErgLPM.f
  ${LIBLOFT}/EM/FromCos/cmigdFunc.f
  ${LIBLOFT}/EM/FromCos/cpairErgLPM.f
  ${LIBLOFT}/EM/FromCos/cphotop.f

  ### LibLoft/EM/FromEpi
  ${LIBLOFT}/EM/FromEpi/epBrCSamp.f
  ${LIBLOFT}/EM/FromEpi/epBrHSamp.f
  ${LIBLOFT}/EM/FromEpi/epBrLSamp.f
  ${LIBLOFT}/EM/FromEpi/epBrSSamp.f
  ${LIBLOFT}/EM/FromEpi/epBrSamp.f
  ${LIBLOFT}/EM/FromEpi/epCompScrBr.f
  ${LIBLOFT}/EM/FromEpi/epCompScrPr.f
  ${LIBLOFT}/EM/FromEpi/epCoulombC.f
  ${LIBLOFT}/EM/FromEpi/epExpot.f
  ${LIBLOFT}/EM/FromEpi/epGetEffZA.f
  ${LIBLOFT}/EM/FromEpi/epKnock.f
  ${LIBLOFT}/EM/FromEpi/epMolBlogB.f90
  ${LIBLOFT}/EM/FromEpi/epMolfunc.f90
  ${LIBLOFT}/EM/FromEpi/epPBA.f
  ${LIBLOFT}/EM/FromEpi/epPrCSamp.f
  ${LIBLOFT}/EM/FromEpi/epPrHSamp.f
  ${LIBLOFT}/EM/FromEpi/epPrLSamp.f
  ${LIBLOFT}/EM/FromEpi/epPrSamp.f
  ${LIBLOFT}/EM/FromEpi/epPrTSamp.f
  ${LIBLOFT}/EM/FromEpi/epReadMTbl.f
  ${LIBLOFT}/EM/FromEpi/epReadTab.f
  ${LIBLOFT}/EM/FromEpi/epReadXXsec.f
  ${LIBLOFT}/EM/FromEpi/epSampMol.f
  ${LIBLOFT}/EM/FromEpi/epSampMolRedA.f90
  ${LIBLOFT}/EM/FromEpi/epSetPhotoE.f
  ${LIBLOFT}/EM/FromEpi/epSetSTblCns.f
  ${LIBLOFT}/EM/FromEpi/epSetTcut.f
  ${LIBLOFT}/EM/FromEpi/epSmpBremAng2BN.f90
  ${LIBLOFT}/EM/FromEpi/epSmpBremAngRelF.f90
  ${LIBLOFT}/EM/FromEpi/epSmpBremAngTsai.f90
  ${LIBLOFT}/EM/FromEpi/epSrim.f90
  ${LIBLOFT}/EM/FromEpi/epStern.f
  ${LIBLOFT}/EM/FromEpi/epUrban.f
  ${LIBLOFT}/EM/FromEpi/epX0.f
  ${LIBLOFT}/EM/FromEpi/epX0Old.f
  ${LIBLOFT}/EM/FromEpi/epXrayp.f
  ${LIBLOFT}/EM/FromEpi/epanihi.f
  ${LIBLOFT}/EM/FromEpi/epbhabha.f
  ${LIBLOFT}/EM/FromEpi/epcompton.f
  ${LIBLOFT}/EM/FromEpi/epdEdxDenC.f
  ${LIBLOFT}/EM/FromEpi/epdedxNone.f
  ${LIBLOFT}/EM/FromEpi/epdedxe.f
  ${LIBLOFT}/EM/FromEpi/epdedxflhv.f
  ${LIBLOFT}/EM/FromEpi/epdedxhvy.f
  ${LIBLOFT}/EM/FromEpi/epdoEMInte.f
  ${LIBLOFT}/EM/FromEpi/epinteElec.f
  ${LIBLOFT}/EM/FromEpi/epinteG.f
  ${LIBLOFT}/EM/FromEpi/epmoller.f
  ${LIBLOFT}/EM/FromEpi/epmulScat.f
  ${LIBLOFT}/EM/FromEpi/epphotoE.f
  ${LIBLOFT}/EM/FromEpi/eprdSmpTbl.f
  ${LIBLOFT}/EM/FromEpi/epsync.f
  ${LIBLOFT}/EM/FromEpi/epwtStern.f

  ### LibLoft/General/Clib
  ${LIBLOFT}/General/Clib/c2bdcy.f
  ${LIBLOFT}/General/Clib/cadjm.f
  ${LIBLOFT}/General/Clib/cbst0.f
  ${LIBLOFT}/General/Clib/cbst1.f
  ${LIBLOFT}/General/Clib/ce2pp.f
  ${LIBLOFT}/General/Clib/cgeqm.f
  ${LIBLOFT}/General/Clib/cgetBeta.f
  ${LIBLOFT}/General/Clib/cgetRotMat3.f
  ${LIBLOFT}/General/Clib/cgetRotMat4.f
  ${LIBLOFT}/General/Clib/cgetcm.f
  ${LIBLOFT}/General/Clib/cgetctau.f
  ${LIBLOFT}/General/Clib/cgetlf.f
  ${LIBLOFT}/General/Clib/cibst1.f
  ${LIBLOFT}/General/Clib/cloreb.f
  ${LIBLOFT}/General/Clib/clorep.f
  ${LIBLOFT}/General/Clib/clorez.f90
  ${LIBLOFT}/General/Clib/cmkptc.f
  ${LIBLOFT}/General/Clib/cnbdcy.f
  ${LIBLOFT}/General/Clib/cpm2e.f
  ${LIBLOFT}/General/Clib/cpxyzp.f
  ${LIBLOFT}/General/Clib/crot3mom.f
  ${LIBLOFT}/General/Clib/crot3vec.f
  ${LIBLOFT}/General/Clib/cscalerProd.f
  ${LIBLOFT}/General/Clib/ctestOnShell.f
  ${LIBLOFT}/General/Clib/ctransVectZ.f
  ${LIBLOFT}/General/Clib/cvecProd.f
  ${LIBLOFT}/General/Clib/cyeta.f

  ### LibLoft/General/FromCos
  ${LIBLOFT}/General/FromCos/csampAF.f90
  ${LIBLOFT}/General/FromCos/k16pGaussLeg.f
  ${LIBLOFT}/General/FromCos/k4ptdi.f
  ${LIBLOFT}/General/FromCos/kNewtonRaphson.f90
  ${LIBLOFT}/General/FromCos/kSNewtonRaphson.f90
  ${LIBLOFT}/General/FromCos/kSRunge-Kutta4.f90
  ${LIBLOFT}/General/FromCos/kSbinChopWP.f90
  ${LIBLOFT}/General/FromCos/kScsplCoef.f
  ${LIBLOFT}/General/FromCos/kScsplDif.f
  ${LIBLOFT}/General/FromCos/kScsplInteg.f
  ${LIBLOFT}/General/FromCos/kScsplIntp.f
  ${LIBLOFT}/General/FromCos/kbetar.f
  ${LIBLOFT}/General/FromCos/kbinChop.f
  ${LIBLOFT}/General/FromCos/kbinChopWP.f90
  ${LIBLOFT}/General/FromCos/kbinom.f
  ${LIBLOFT}/General/FromCos/kceles.f
  ${LIBLOFT}/General/FromCos/kcombi.f
  ${LIBLOFT}/General/FromCos/kcossn.f
  ${LIBLOFT}/General/FromCos/kcountFields.f90
  ${LIBLOFT}/General/FromCos/kcsplCoef.f
  ${LIBLOFT}/General/FromCos/kcsplDif.f
  ${LIBLOFT}/General/FromCos/kcsplInteg.f
  ${LIBLOFT}/General/FromCos/kcsplIntp.f
  ${LIBLOFT}/General/FromCos/kcsr1idx.f
  ${LIBLOFT}/General/FromCos/kcubicEq.f90
  ${LIBLOFT}/General/FromCos/kdexpIntF.f
  ${LIBLOFT}/General/FromCos/kdexpIntFWP.f90
  ${LIBLOFT}/General/FromCos/kdexpIntFb.f
  ${LIBLOFT}/General/FromCos/kdir2deg.f
  ${LIBLOFT}/General/FromCos/kdmachmnmx.f
  ${LIBLOFT}/General/FromCos/kdwhereis.f
  ${LIBLOFT}/General/FromCos/kerf.f
  ${LIBLOFT}/General/FromCos/kerfc.f
  ${LIBLOFT}/General/FromCos/kexpiC.f
  ${LIBLOFT}/General/FromCos/kexpot.f
  ${LIBLOFT}/General/FromCos/kfrge.f
  ${LIBLOFT}/General/FromCos/kgamma.f
  ${LIBLOFT}/General/FromCos/kgauss.f
  ${LIBLOFT}/General/FromCos/kggamma.f90
  ${LIBLOFT}/General/FromCos/kgsub.f
  ${LIBLOFT}/General/FromCos/kintp3.f
  ${LIBLOFT}/General/FromCos/klena.f
  ${LIBLOFT}/General/FromCos/kmanbit.f
  ${LIBLOFT}/General/FromCos/knbino.f
  ${LIBLOFT}/General/FromCos/kpmnx.f
  ${LIBLOFT}/General/FromCos/kpoisn.f
  ${LIBLOFT}/General/FromCos/kpolintp.f
  ${LIBLOFT}/General/FromCos/kpolintp2.f
  ${LIBLOFT}/General/FromCos/kpolintp2S.f
  ${LIBLOFT}/General/FromCos/kpolintpS.f
  ${LIBLOFT}/General/FromCos/kpolygC.f
  ${LIBLOFT}/General/FromCos/kqsort.f
  ${LIBLOFT}/General/FromCos/kqsortc.f
  ${LIBLOFT}/General/FromCos/kqsortd.f
  ${LIBLOFT}/General/FromCos/kqsorth.f
  ${LIBLOFT}/General/FromCos/kqsorti.f
  ${LIBLOFT}/General/FromCos/kqsortr.f
  ${LIBLOFT}/General/FromCos/kquadEq.f90
  ${LIBLOFT}/General/FromCos/kroundrv.f
  ${LIBLOFT}/General/FromCos/ksampAF.f
  ${LIBLOFT}/General/FromCos/ksampForBPA.f90
  ${LIBLOFT}/General/FromCos/ksampLin.f
  ${LIBLOFT}/General/FromCos/ksampPEang.f
  ${LIBLOFT}/General/FromCos/ksampPw.f
  ${LIBLOFT}/General/FromCos/ksampRSA.f
  ${LIBLOFT}/General/FromCos/ksbwig.f
  ${LIBLOFT}/General/FromCos/kseblk.f
  ${LIBLOFT}/General/FromCos/ksgamd.f
  ${LIBLOFT}/General/FromCos/ksmooth.f
  ${LIBLOFT}/General/FromCos/ksmpConExp.f
  ${LIBLOFT}/General/FromCos/ksortinv.f
  ${LIBLOFT}/General/FromCos/ksplandau.f
  ${LIBLOFT}/General/FromCos/ksplit.f
  ${LIBLOFT}/General/FromCos/ksupblank.f
  ${LIBLOFT}/General/FromCos/ktrpzIntT2.f
  ${LIBLOFT}/General/FromCos/ktrpzIntegT.f
  ${LIBLOFT}/General/FromCos/kvec_prod3.f90
  ${LIBLOFT}/General/FromCos/kzeta2a.f
  ${LIBLOFT}/General/FromCos/rnd.f

  ### LibLoft/General/FromEpi
  ${LIBLOFT}/General/FromEpi/afsep.f
  ${LIBLOFT}/General/FromEpi/epcos2dir.f90
  ${LIBLOFT}/General/FromEpi/episoAngle.f
  ${LIBLOFT}/General/FromEpi/k16pGausLeg2.f
  ${LIBLOFT}/General/FromEpi/kCheckCoPlan.f90
  ${LIBLOFT}/General/FromEpi/kProjAReaOfTriang.f90
  ${LIBLOFT}/General/FromEpi/kbchop.f
  ${LIBLOFT}/General/FromEpi/kbivarNorm.f
  ${LIBLOFT}/General/FromEpi/kbtest.f
  ${LIBLOFT}/General/FromEpi/kcquadEq.f90
  ${LIBLOFT}/General/FromEpi/kcsr1.f
  ${LIBLOFT}/General/FromEpi/kdgtest.f
  ${LIBLOFT}/General/FromEpi/kexpi.f
  ${LIBLOFT}/General/FromEpi/kgetAngleBtw2Planes.f90
  ${LIBLOFT}/General/FromEpi/kgetNormalVec.f90
  ${LIBLOFT}/General/FromEpi/kinout.f90
  ${LIBLOFT}/General/FromEpi/kinout2.f90
  ${LIBLOFT}/General/FromEpi/kinout3.f90
  ${LIBLOFT}/General/FromEpi/kisInPolyhed0.f90
  ${LIBLOFT}/General/FromEpi/kklib.f
  ${LIBLOFT}/General/FromEpi/kphiinside.f90
  ${LIBLOFT}/General/FromEpi/kpolyg.f
  ${LIBLOFT}/General/FromEpi/kquarticEq.f90
  ${LIBLOFT}/General/FromEpi/krmLeadingB.f90
  ${LIBLOFT}/General/FromEpi/ksbeta.f
  ${LIBLOFT}/General/FromEpi/kseeinbox.f
  ${LIBLOFT}/General/FromEpi/kseeincyl.f
  ${LIBLOFT}/General/FromEpi/ksmpintbeta.f
  ${LIBLOFT}/General/FromEpi/ksort1.f
  ${LIBLOFT}/General/FromEpi/ksplit2.f90
  ${LIBLOFT}/General/FromEpi/ksx2.f
  ${LIBLOFT}/General/FromEpi/kxpLand3vpB.f90
  ${LIBLOFT}/General/FromEpi/kxpLineAndTorus.f90
  ${LIBLOFT}/General/FromEpi/kxplCylWall.f90
  ${LIBLOFT}/General/FromEpi/kxplPipe.f90
  ${LIBLOFT}/General/FromEpi/kxplPolyhed0.f90
  ${LIBLOFT}/General/FromEpi/kxplccl.f
  ${LIBLOFT}/General/FromEpi/kxplhorse.f
  ${LIBLOFT}/General/FromEpi/kxplsq.f
  ${LIBLOFT}/General/FromEpi/kxplvsq.f90
  ${LIBLOFT}/General/FromEpi/psiim.f

  ### LibLoft/Had/Event
  ${LIBLOFT}/Had/Event/Aacol/caveInteNuc.f
  ${LIBLOFT}/Had/Event/Aacol/cheavyInt.f
  ${LIBLOFT}/Had/Event/GammaN/cgpHad.f
  ${LIBLOFT}/Had/Event/GammaN/csPiAngOfPiN.f
  ${LIBLOFT}/Had/Event/Hacol/caveNoOfSucC.f
  ${LIBLOFT}/Had/Event/Hacol/cfxTgtChg.f
  ${LIBLOFT}/Had/Event/Hacol/chAcol.f
  ${LIBLOFT}/Had/Event/Hacol/chAcolAdhoc.f
  ${LIBLOFT}/Had/Event/Hacol/csampCollInA.f
  ${LIBLOFT}/Had/Event/Hacol/csampFermiM.f
  ${LIBLOFT}/Had/Event/Hacol/csmpColInA2.f
  ${LIBLOFT}/Had/Event/Hncol/cconsvChg.f
  ${LIBLOFT}/Had/Event/Hncol/cgnlp.f
  ${LIBLOFT}/Had/Event/Hncol/cgrap.f
  ${LIBLOFT}/Had/Event/Hncol/chncol.f
  ${LIBLOFT}/Had/Event/Hncol/cinipipx.f
  ${LIBLOFT}/Had/Event/Hncol/cinippx.f
  ${LIBLOFT}/Had/Event/Hncol/cs2lp.f
  ${LIBLOFT}/Had/Event/Hncol/csampPolAng.f
  ${LIBLOFT}/Had/Event/Hncol/csnchp.f
  ${LIBLOFT}/Had/Event/Hncol/csptxy.f
  ${LIBLOFT}/Had/Event/Hncol/cspwpt.f
  ${LIBLOFT}/Had/Event/cNucMass.f90
  ${LIBLOFT}/Had/Event/cTargetFrag.f90

  ### LibLoft/Had/Import/DPM
  ${LIBLOFT}/Had/Import/DPM/Interface/cdpmjet.f
  ${LIBLOFT}/Had/Import/DPM/dpmjet3.f
  ${LIBLOFT}/Had/Import/DPM/phojet.f
  ${LIBLOFT}/Had/Import/DPM/pythia6115.f

  ### LibLoft/Had/Import/EOPS
  ${LIBLOFT}/Had/Import/EPOS/Interface/cepos.f
  ${LIBLOFT}/Had/Import/EPOS/Interface/cepos2cos.f90
  ${LIBLOFT}/Had/Import/EPOS/Interface/cepos2cosB.f
  ${LIBLOFT}/Had/Import/EPOS/epos-bas-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-con-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-dky-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-dro-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-ems-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-fra-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-hnb-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-ico-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-ids-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-int-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-jps-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-omg-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-par-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-qsh-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-rsh-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-sem-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-tim-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-uti-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-xan-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos-xpr-lhc-v3700.f
  ${LIBLOFT}/Had/Import/EPOS/epos_j.f
  ${LIBLOFT}/Had/Import/EPOS/eposm.f
  ${LIBLOFT}/Had/Import/EPOS/eposu_no.f

  ### LibLoft/Had/Import/Fritiof
  ${LIBLOFT}/Had/Import/Fritiof/Interface/cfrevent.f
  ${LIBLOFT}/Had/Import/Fritiof/Interface/chALund.f
  ${LIBLOFT}/Had/Import/Fritiof/aadn.f
  ${LIBLOFT}/Had/Import/Fritiof/aado.f
  ${LIBLOFT}/Had/Import/Fritiof/aath.f
  ${LIBLOFT}/Had/Import/Fritiof/aatj.f

  ### LibLoft/Had/Import/Gheisha
  ${LIBLOFT}/Had/Import/Gheisha/Interface/ccos2gheCode.f
  ${LIBLOFT}/Had/Import/Gheisha/Interface/cghe2cos.f
  ${LIBLOFT}/Had/Import/Gheisha/Interface/cgheDecay.f
  ${LIBLOFT}/Had/Import/Gheisha/Interface/chAGheisha.f
  ${LIBLOFT}/Had/Import/Gheisha/Interface/cxsecGheisha.f
  ${LIBLOFT}/Had/Import/Gheisha/acfsin.f
  ${LIBLOFT}/Had/Import/Gheisha/add.f
  ${LIBLOFT}/Had/Import/Gheisha/ang.f
  ${LIBLOFT}/Had/Import/Gheisha/atomas.f
  ${LIBLOFT}/Had/Import/Gheisha/captur.f
  ${LIBLOFT}/Had/Import/Gheisha/casal0.f
  ${LIBLOFT}/Had/Import/Gheisha/casaom.f
  ${LIBLOFT}/Had/Import/Gheisha/casasm.f
  ${LIBLOFT}/Had/Import/Gheisha/casasp.f
  ${LIBLOFT}/Had/Import/Gheisha/casax0.f
  ${LIBLOFT}/Had/Import/Gheisha/casaxm.f
  ${LIBLOFT}/Had/Import/Gheisha/casfrg.f
  ${LIBLOFT}/Had/Import/Gheisha/cask0.f
  ${LIBLOFT}/Had/Import/Gheisha/cask0b.f
  ${LIBLOFT}/Had/Import/Gheisha/caskm.f
  ${LIBLOFT}/Had/Import/Gheisha/caskp.f
  ${LIBLOFT}/Had/Import/Gheisha/casl0.f
  ${LIBLOFT}/Had/Import/Gheisha/casn.f
  ${LIBLOFT}/Had/Import/Gheisha/casnb.f
  ${LIBLOFT}/Had/Import/Gheisha/casom.f
  ${LIBLOFT}/Had/Import/Gheisha/casp.f
  ${LIBLOFT}/Had/Import/Gheisha/caspb.f
  ${LIBLOFT}/Had/Import/Gheisha/caspim.f
  ${LIBLOFT}/Had/Import/Gheisha/caspip.f
  ${LIBLOFT}/Had/Import/Gheisha/cassm.f
  ${LIBLOFT}/Had/Import/Gheisha/cassp.f
  ${LIBLOFT}/Had/Import/Gheisha/casx0.f
  ${LIBLOFT}/Had/Import/Gheisha/casxm.f
  ${LIBLOFT}/Had/Import/Gheisha/chanwx.f
  ${LIBLOFT}/Had/Import/Gheisha/cinema.f
  ${LIBLOFT}/Had/Import/Gheisha/compo.f
  ${LIBLOFT}/Had/Import/Gheisha/coranh.f
  ${LIBLOFT}/Had/Import/Gheisha/coscat.f
  ${LIBLOFT}/Had/Import/Gheisha/defs.f
  ${LIBLOFT}/Had/Import/Gheisha/defs1.f
  ${LIBLOFT}/Had/Import/Gheisha/dlpnuc.f
  ${LIBLOFT}/Had/Import/Gheisha/dlpsor.f
  ${LIBLOFT}/Had/Import/Gheisha/dotes2.f
  ${LIBLOFT}/Had/Import/Gheisha/dotnuc.f
  ${LIBLOFT}/Had/Import/Gheisha/dpdnuc.f
  ${LIBLOFT}/Had/Import/Gheisha/exnu.f
  ${LIBLOFT}/Had/Import/Gheisha/fctcos.f
  ${LIBLOFT}/Had/Import/Gheisha/fermi.f
  ${LIBLOFT}/Had/Import/Gheisha/fissio.f
  ${LIBLOFT}/Had/Import/Gheisha/gdpdk.f
  ${LIBLOFT}/Had/Import/Gheisha/genxpt.f
  ${LIBLOFT}/Had/Import/Gheisha/gheini.f
  ${LIBLOFT}/Had/Import/Gheisha/gheish.f
  ${LIBLOFT}/Had/Import/Gheisha/ghstop.f
  ${LIBLOFT}/Had/Import/Gheisha/ghtype.f
  ${LIBLOFT}/Had/Import/Gheisha/gmate.f
  ${LIBLOFT}/Had/Import/Gheisha/gnslwd.f
  ${LIBLOFT}/Had/Import/Gheisha/gpart.f
  ${LIBLOFT}/Had/Import/Gheisha/gpdk.f
  ${LIBLOFT}/Had/Import/Gheisha/gpghei.f
  ${LIBLOFT}/Had/Import/Gheisha/gpions.f
  ${LIBLOFT}/Had/Import/Gheisha/gspart.f
  ${LIBLOFT}/Had/Import/Gheisha/gutil.f
  ${LIBLOFT}/Had/Import/Gheisha/hadini.f
  ${LIBLOFT}/Had/Import/Gheisha/hadnuc.f
  ${LIBLOFT}/Had/Import/Gheisha/kmabs.f
  ${LIBLOFT}/Had/Import/Gheisha/lengtx.f
  ${LIBLOFT}/Had/Import/Gheisha/lor.f
  ${LIBLOFT}/Had/Import/Gheisha/minv.f
  ${LIBLOFT}/Had/Import/Gheisha/nbanh.f
  ${LIBLOFT}/Had/Import/Gheisha/nfac.f
  ${LIBLOFT}/Had/Import/Gheisha/normal.f
  ${LIBLOFT}/Had/Import/Gheisha/nucrec.f
  ${LIBLOFT}/Had/Import/Gheisha/pbanh.f
  ${LIBLOFT}/Had/Import/Gheisha/pcsdat.f
  ${LIBLOFT}/Had/Import/Gheisha/pdgtyp.f
  ${LIBLOFT}/Had/Import/Gheisha/phasp.f
  ${LIBLOFT}/Had/Import/Gheisha/phpnuc.f
  ${LIBLOFT}/Had/Import/Gheisha/pimabs.f
  ${LIBLOFT}/Had/Import/Gheisha/pmltpc.f
  ${LIBLOFT}/Had/Import/Gheisha/poisso.f
  ${LIBLOFT}/Had/Import/Gheisha/ranres.f
  ${LIBLOFT}/Had/Import/Gheisha/rotes2.f
  ${LIBLOFT}/Had/Import/Gheisha/rtmi.f
  ${LIBLOFT}/Had/Import/Gheisha/selfab.f
  ${LIBLOFT}/Had/Import/Gheisha/setcur.f
  ${LIBLOFT}/Had/Import/Gheisha/settrk.f
  ${LIBLOFT}/Had/Import/Gheisha/sginit.f
  ${LIBLOFT}/Had/Import/Gheisha/steep.f
  ${LIBLOFT}/Had/Import/Gheisha/steeq.f
  ${LIBLOFT}/Had/Import/Gheisha/steer.f
  ${LIBLOFT}/Had/Import/Gheisha/stees.f
  ${LIBLOFT}/Had/Import/Gheisha/stpair.f
  ${LIBLOFT}/Had/Import/Gheisha/sub.f
  ${LIBLOFT}/Had/Import/Gheisha/tdelay.f
  ${LIBLOFT}/Had/Import/Gheisha/trac.f
  ${LIBLOFT}/Had/Import/Gheisha/twob.f
  ${LIBLOFT}/Had/Import/Gheisha/twoclu.f

  ### LibLoft/Had/Import/InclusiveDPM
  ${LIBLOFT}/Had/Import/InclusiveDPM/Interface/cincdpm3.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-common/int-crosssection.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-common/int-engine.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-common/int-initvarip.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-distbndry.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-multi.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-sample-cospt.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-v1.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-v2.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-v3.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-v4.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-v5.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-v6.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-v7.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-v8.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-v9.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-variblock.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-vinterface2.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-x1.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-x2.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-x3.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-x4.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-x5.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-x6.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-x7.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-x8.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-x9.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/int-dpmjet3/int-xinterface2.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/math/math-bspl.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/math/math-character.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/math/math-heapsort.f
  ${LIBLOFT}/Had/Import/InclusiveDPM/math/math-lorentz.f

  ### LibLoft/Had/Import/Jam
  ${LIBLOFT}/Had/Import/Jam/Interface/cfragment.f
  ${LIBLOFT}/Had/Import/Jam/Interface/cjam.f
  ${LIBLOFT}/Had/Import/Jam/gfortran.f
  ${LIBLOFT}/Had/Import/Jam/jam.f
  ${LIBLOFT}/Had/Import/Jam/jamana.f
  ${LIBLOFT}/Had/Import/Jam/jambuu.f
  ${LIBLOFT}/Had/Import/Jam/jamcoll.f
  ${LIBLOFT}/Had/Import/Jam/jamcross.f
  ${LIBLOFT}/Had/Import/Jam/jamdat.f
  ${LIBLOFT}/Had/Import/Jam/jamdec.f
  ${LIBLOFT}/Had/Import/Jam/jamdummy.f
  ${LIBLOFT}/Had/Import/Jam/jamhard.f
  ${LIBLOFT}/Had/Import/Jam/jamhij.f
  ${LIBLOFT}/Had/Import/Jam/jampdf.f
  ${LIBLOFT}/Had/Import/Jam/jampyjet.f
  ${LIBLOFT}/Had/Import/Jam/jampysigh.f
  ${LIBLOFT}/Had/Import/Jam/jampythia.f
  ${LIBLOFT}/Had/Import/Jam/jamsoft.f

  ### LibLoft/Had/Import/NUCFRG
  ${LIBLOFT}/Had/Import/NUCFRG/Interface/csampFermiMom.f90
  ${LIBLOFT}/Had/Import/NUCFRG/Interface/csampHFMom.f90
  ${LIBLOFT}/Had/Import/NUCFRG/Interface/csampNucFrag.f
  ${LIBLOFT}/Had/Import/NUCFRG/abbg.f

  ### LibLoft/Had/Import/PHITS
  ${LIBLOFT}/Had/Import/PHITS/Bertini/cbert.f
  ${LIBLOFT}/Had/Import/PHITS/Com/kftp.f
  ${LIBLOFT}/Had/Import/PHITS/Com/ndata01.f
  ${LIBLOFT}/Had/Import/PHITS/Com/phitsrn.f
  ${LIBLOFT}/Had/Import/PHITS/Com/utl01.f
  ${LIBLOFT}/Had/Import/PHITS/Evapo/cnevap.f
  ${LIBLOFT}/Had/Import/PHITS/Interface/cbertini.f
  ${LIBLOFT}/Had/Import/PHITS/Interface/cjqmd.f
  ${LIBLOFT}/Had/Import/PHITS/Interface/cnelas.f
  ${LIBLOFT}/Had/Import/PHITS/Interface/cphits.f
  ${LIBLOFT}/Had/Import/PHITS/JQMD/jqmd.f
  ${LIBLOFT}/Had/Import/PHITS/Nela/nelst.f
  ${LIBLOFT}/Had/Import/PHITS/Xsec/phitsXsec.f
  ${LIBLOFT}/Had/Import/PHITS/bert-bl0.f
  ${LIBLOFT}/Had/Import/PHITS/bert-bl1.f
  ${LIBLOFT}/Had/Import/PHITS/bert-bl2.f
  ${LIBLOFT}/Had/Import/PHITS/cphits2cos.f
  ${LIBLOFT}/Had/Import/PHITS/cphitsOut.f
  ${LIBLOFT}/Had/Import/PHITS/cprePhits.f
  ${LIBLOFT}/Had/Import/PHITS/gemset.f
  ${LIBLOFT}/Had/Import/PHITS/ndata01blk.f
  ${LIBLOFT}/Had/Import/PHITS/qmddflt.f

  ### LibLoft/Had/Import/QGS
  ${LIBLOFT}/Had/Import/QGS/Interface/cQGSjet.f
  ${LIBLOFT}/Had/Import/QGS/qgprox.f
  ${LIBLOFT}/Had/Import/QGS/qgsjet-II-04.f

  ### LibLoft/Had/Import/Sibyll
  ${LIBLOFT}/Had/Import/Sibyll/GASDEV.f
  ${LIBLOFT}/Had/Import/Sibyll/Interface/csibSetStblPtcl.f
  ${LIBLOFT}/Had/Import/Sibyll/Interface/csibyll.f
  ${LIBLOFT}/Had/Import/Sibyll/Interface/csibyllXs.f
  ${LIBLOFT}/Had/Import/Sibyll/Interface/csibyllXsAAir.f
  ${LIBLOFT}/Had/Import/Sibyll/sibyll2.3c02.f
  ${LIBLOFT}/Had/Import/Sofia/Interface/csofia.f
  ${LIBLOFT}/Had/Import/Sofia/jetset74dp.f
  ${LIBLOFT}/Had/Import/Sofia/rndm.f
  ${LIBLOFT}/Had/Import/Sofia/sofia.f

  ### LibLoft/Had/Interface
  ${LIBLOFT}/Had/Interface/cfixModel.f
  ${LIBLOFT}/Had/Interface/cintModels.f
  ${LIBLOFT}/Had/Interface/cinteNEP.f
  ${LIBLOFT}/Had/Interface/cinteNEPsub.f

  ### LibLoft/Had/Xsection
  ${LIBLOFT}/Had/Xsection/cAAXsec2.f
  ${LIBLOFT}/Had/Xsection/cGetXsec.f
  ${LIBLOFT}/Had/Xsection/cgpXsec.f
  ${LIBLOFT}/Had/Xsection/cgpxs1.f
  ${LIBLOFT}/Had/Xsection/cinelx.f
  ${LIBLOFT}/Had/Xsection/ckmnTotXs.f
  ${LIBLOFT}/Had/Xsection/ckmpTotXs.f
  ${LIBLOFT}/Had/Xsection/ckpnTotXs.f
  ${LIBLOFT}/Had/Xsection/ckppTotXs.f
  ${LIBLOFT}/Had/Xsection/cnbarpTotXs.f
  ${LIBLOFT}/Had/Xsection/cnpTotXs.f
  ${LIBLOFT}/Had/Xsection/cpbarpTotXs.f
  ${LIBLOFT}/Had/Xsection/cpdgXs.f
  ${LIBLOFT}/Had/Xsection/cpdgXs200GeV.f90
  ${LIBLOFT}/Had/Xsection/cpimnTotXs.f
  ${LIBLOFT}/Had/Xsection/cpimpTotXs.f
  ${LIBLOFT}/Had/Xsection/cpipnTotXs.f
  ${LIBLOFT}/Had/Xsection/cpippTotXs.f
  ${LIBLOFT}/Had/Xsection/cppTotXs.f
  ${LIBLOFT}/Had/Xsection/cppXsec.f
  ${LIBLOFT}/Had/Xsection/cprotonAXsec.f
  ${LIBLOFT}/Had/Xsection/ctotx.f
  ${LIBLOFT}/Had/Xsection/cxAbyxpXsec.f
  ${LIBLOFT}/Had/Xsection/cxp2xAXsec.f
  ${LIBLOFT}/Had/Xsection/cxsSpecial.f

  ### LibLoft/Manager
  ${LIBLOFT}/Manager/Ccode/kgetenv.c
  ${LIBLOFT}/Manager/c2lowerCase.f
  ${LIBLOFT}/Manager/cCodeConv.f
  ${LIBLOFT}/Manager/cepControl.f
  ${LIBLOFT}/Manager/cerrorMsg.f
  ${LIBLOFT}/Manager/cgetLoginN.f
  ${LIBLOFT}/Manager/ckf2cos.f
  ${LIBLOFT}/Manager/cmkSeed.f
  ${LIBLOFT}/Manager/copenf.f
  ${LIBLOFT}/Manager/creadParam.f
  ${LIBLOFT}/Manager/cskipComment.f
  ${LIBLOFT}/Manager/cwriteParam.f
  ${LIBLOFT}/Manager/cwriteSeed.f

  ### LibLoft/Mu
  ${LIBLOFT}/Mu/FromEpi/epRdmuTab.f
  ${LIBLOFT}/Mu/FromEpi/epSetMu.f
  ${LIBLOFT}/Mu/FromEpi/epSetmuSTab.f
  ${LIBLOFT}/Mu/FromEpi/epmuAuxP.f
  ${LIBLOFT}/Mu/FromEpi/epmuBrsmp.f
  ${LIBLOFT}/Mu/FromEpi/epmuNsmp.f
  ${LIBLOFT}/Mu/FromEpi/epmuPrsmp.f
  ${LIBLOFT}/Mu/FromEpi/epmuSetCnst.f
  ${LIBLOFT}/Mu/FromEpi/epmudEdx.f
  ${LIBLOFT}/Mu/FromEpi/epmudsdv.f
  ${LIBLOFT}/Mu/FromEpi/epmudsdvdr.f
  ${LIBLOFT}/Mu/FromEpi/epmuvmax.f
  ${LIBLOFT}/Mu/FromEpi/epmuvmax2.f
  ${LIBLOFT}/Mu/FromEpi/eprdmuSmpTbl.f
  ${LIBLOFT}/Mu/cmuInte.f

  ### LibLoft/PreIntre
  ${LIBLOFT}/PreInte/cdecayLeng.f
  ${LIBLOFT}/PreInte/cepSampEIntL.f
  ${LIBLOFT}/PreInte/cepSampGIntL.f
  ${LIBLOFT}/PreInte/cepSampNEPIntL.f
  ${LIBLOFT}/PreInte/cseeColPossible.f
  ${LIBLOFT}/PreInte/epprocAlias.f
  ${LIBLOFT}/PreInte/eprdMFile.f

  ### LibLoft/Primary
  ${LIBLOFT}/Primary/ciniSPrim.f
  ${LIBLOFT}/Primary/cmkPrimSTbl.f
  ${LIBLOFT}/Primary/cprintPrim.f
  ${LIBLOFT}/Primary/cprocPrimDt.f
  ${LIBLOFT}/Primary/crdPrimData.f
  ${LIBLOFT}/Primary/crigCut.f
  ${LIBLOFT}/Primary/csampPrimary.f


  ### Cosmos/Atmosphere
  ${COSMOSTOP}/Atmosphere/AtmosModel/ciniAtmos.f
  ${COSMOSTOP}/Atmosphere/AtmosModel/creadAtmosD.f
  ${COSMOSTOP}/Atmosphere/AtmosModel/cstdatmos0.f
  ${COSMOSTOP}/Atmosphere/catmosCnst1.f
  ${COSMOSTOP}/Atmosphere/catmosDef.f90
  ${COSMOSTOP}/Atmosphere/catmosrho.f
  ${COSMOSTOP}/Atmosphere/catmosutil.f
  ${COSMOSTOP}/Atmosphere/ceCent2llh.f
  ${COSMOSTOP}/Atmosphere/cgetMedia.f
  ${COSMOSTOP}/Atmosphere/cl2tTbl.f
  ${COSMOSTOP}/Atmosphere/cl2tTblA.f
  ${COSMOSTOP}/Atmosphere/clen2thick.f
  ${COSMOSTOP}/Atmosphere/cllh2eCent.f
  ${COSMOSTOP}/Atmosphere/cllh2sph.f
  ${COSMOSTOP}/Atmosphere/csph2llh.f
  ${COSMOSTOP}/Atmosphere/cthick2len.f

  ### Cosmos/Manager
  ${COSMOSTOP}/Manager/cbeginRun.f
  ${COSMOSTOP}/Manager/cendRun.f
  ${COSMOSTOP}/Manager/ceventLoop.f
  ${COSMOSTOP}/Manager/cmanager.f
  ${COSMOSTOP}/Manager/cqEventNo.f
  ${COSMOSTOP}/Manager/cqUHookr.f
  ${COSMOSTOP}/Manager/cqverstion.f

  ### Cosmos/Tracking
  ${COSMOSTOP}/Tracking/AS/cNeByApproxB.f
  ${COSMOSTOP}/Tracking/AS/cavedEdx.f
  ${COSMOSTOP}/Tracking/AS/cgetMoliereU.f
  ${COSMOSTOP}/Tracking/AS/cgetNmu.f
  ${COSMOSTOP}/Tracking/AS/clatDist.f
  ${COSMOSTOP}/Tracking/AS/cobAS.f
  ${COSMOSTOP}/Tracking/Geomag/cbDefByRK.f
  ${COSMOSTOP}/Tracking/Geomag/cgeomag.f
  ${COSMOSTOP}/Tracking/Geomag/cmagneticDef.f
  ${COSMOSTOP}/Tracking/Geomag/cprintMagF.f
  ${COSMOSTOP}/Tracking/Geomag/csetMagField.f
  ${COSMOSTOP}/Tracking/Geomag/ctransMagTo.f
  ${COSMOSTOP}/Tracking/c3DV2DDCos.f
  ${COSMOSTOP}/Tracking/c4to3mom.f
  ${COSMOSTOP}/Tracking/cEfield0.f
  ${COSMOSTOP}/Tracking/ccompPathEnd.f
  ${COSMOSTOP}/Tracking/cdefByMagAndE.f
  ${COSMOSTOP}/Tracking/cdet2xyzFE.f
  ${COSMOSTOP}/Tracking/ce2p.f
  ${COSMOSTOP}/Tracking/ceCent2sph.f
  ${COSMOSTOP}/Tracking/cexcessLen.f
  ${COSMOSTOP}/Tracking/cfinTracking.f
  ${COSMOSTOP}/Tracking/cfixProc.f
  ${COSMOSTOP}/Tracking/cgetMFP.f
  ${COSMOSTOP}/Tracking/cgetZenith.f
  ${COSMOSTOP}/Tracking/cifDead.f
  ${COSMOSTOP}/Tracking/cifXObsSite.f
  ${COSMOSTOP}/Tracking/ciniSmpIntL.f
  ${COSMOSTOP}/Tracking/ciniTracking.f
  ${COSMOSTOP}/Tracking/cinitObs.f
  ${COSMOSTOP}/Tracking/cinteraction.f
  ${COSMOSTOP}/Tracking/cmaxMovLen.f
  ${COSMOSTOP}/Tracking/cmkIncident.f
  ${COSMOSTOP}/Tracking/cobservation.f
  ${COSMOSTOP}/Tracking/cprObsSite.f
  ${COSMOSTOP}/Tracking/cpreSampEIntL.f
  ${COSMOSTOP}/Tracking/cpreSampGIntL.f
  ${COSMOSTOP}/Tracking/cpreSampNEPIntL.f
  ${COSMOSTOP}/Tracking/cprintCoord.f
  ${COSMOSTOP}/Tracking/cprintObs.f
  ${COSMOSTOP}/Tracking/cputCerenkov.f
  ${COSMOSTOP}/Tracking/cputTrInfo.f
  ${COSMOSTOP}/Tracking/cresetDirec.f
  ${COSMOSTOP}/Tracking/cresetMom.f
  ${COSMOSTOP}/Tracking/csPrimAng.f
  ${COSMOSTOP}/Tracking/csetCoord.f
  ${COSMOSTOP}/Tracking/csetPos.f
  ${COSMOSTOP}/Tracking/csph2eCent.f
  ${COSMOSTOP}/Tracking/cstack.f
  ${COSMOSTOP}/Tracking/ctracking.f
  ${COSMOSTOP}/Tracking/ctransCoord2.f
  ${COSMOSTOP}/Tracking/ctransVectZx.f
  ${COSMOSTOP}/Tracking/cvec2dcos.f
  ${COSMOSTOP}/Tracking/cxyz2det.f
  ${COSMOSTOP}/Tracking/cxyz2prim.f
  )

