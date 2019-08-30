#include "BlockData/cblkGene.h"
!  #include "jamdat.f"
#include "jamdummy.h"
#include "phitsdummy.h"
! #include "Zphitsblk.h"
!       main program of cosmos
      program Cosmos
#include "ZcosmosExt.h"
#include "jamextrn.inc"
!///////////
!         call mydummy
!         call myPhitsDummy
!/////////////
         call cmanager   ! call Manager for Cosmos Simulation
      end


