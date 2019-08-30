#include "Zintmodel.h"
      subroutine getImpactParam(bout)
      implicit none
#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"

      real(8),intent(out):: bout  

      real(8):: RPROJ, RTARG,BIMPAC
      integer::       NWTSAM,NWASAM,NWBSAM,NWTACC,NWAACC,NWBACC
      COMMON /DTGLCP/ RPROJ,RTARG,BIMPAC,
     &                NWTSAM,NWASAM,NWBSAM,NWTACC,NWAACC,NWBACC
      integer,save::count=0

!       for qgsjet-II/I
      integer,parameter:: iapmax=209
      real(8):: xa, xb, b
      common /qgarr7/  xa(iapmax,3),xb(iapmax,3),b

      if( ActiveMdl == "dpmjet3" ) then
         bout = BIMPAC
      elseif( ActiveMdl(1:6) == "qgsjet" ) then
         bout = b
      elseif( ActiveMdl == "epos" ) then
         call getEposB(bout) 
      elseif (ActiveMdl == "sibyll" ) then
         call getSibyllB(bout)
      else
         if( count < 10 ) then
            count =count + 1
            write(0,*)
     *      ' IntModel =', ActiveMdl, ' not ready for getImpactParam' 
            write(0,*) 'so  -1.0 is returned'
         endif
         bout = -1.0
      endif
      end
      subroutine getEposB(bout)
      implicit none
      real(8),intent(out):: bout
!---------------------------------------------------------------------------
!                   epos event common block
!---------------------------------------------------------------------------

      real        phievt,bimevt,pmxevt,egyevt
     *,xbjevt,qsqevt,zppevt,zptevt
      integer     nevt,kolevt,koievt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,nglevt,minfra,maxfra
      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra

      bout = bimevt
      end
      subroutine getSibyllB(bout)
      implicit none
      real(8),intent(out):: bout

      integer:: IAMAX, NTRY, NA
      real:: B, BMAX
      PARAMETER (IAMAX=56)
      COMMON /S_CNCM0/ B, BMAX, NTRY, NA
      bout = B
      end
