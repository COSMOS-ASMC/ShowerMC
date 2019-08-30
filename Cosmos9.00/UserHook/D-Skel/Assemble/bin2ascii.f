      implicit none
#define  FNODATDEF 33
!          echo xxx.bdat | bin2ascii...
!
#include "../FleshHist/Zprivate1.h"
#include "../FleshHist/Zprivate3.h"
      character*120 file
      real*8  Et, wx, wy, wz, firstdep
      integer EventNo
      integer*2  code
      integer i

      read(*,'(a)') file
      write(0,*) ' file name is '
      write(0,*) file
      open(fnodat, file=file, form="unformatted")
!     read(fnodat, err=11)
!     *        EventNo, code,   Et,  wx, wy, wz, firstdep
!      goto 12
! 11   continue
!      write(0,*) ' assuming old data without firstdep info'
!       old version data has no firstdep
!      rewind(fnodat)
!      read(fnodat)
!     *        EventNo, code,   Et,  wx, wy, wz
!      firstdep=-1. 
! 12   continue
!      write(*,'("i ",  i3,  i4, 1pE12.3,0p3f11.7, f8.1)')
!     *        EventNo, code,   Et,  wx, wy, wz, firstdep
      
      do while(.true.)
         read(fnodat,end=100) bufc, (buf(i), i=1, bufc)
         do i = 1, bufc
#if KeepWeight != yes
            write(*,
     *      '(6i3, 1pE11.3, 0p,f6.1,1p2E11.3,0p, 2f8.4,f10.6)')
     *        buf(i).ldep,  buf(i).code,  buf(i).subcode,
     *        buf(i).charge, buf(i).ridx, buf(i).faiidx,
     *        buf(i).rinmu, buf(i).fai, buf(i).Ek,
     *        buf(i).t, buf(i).wx, buf(i).wy, buf(i).wz
#else
            write(*,
     *   '(6i3, 1pE11.3, 0p,f6.1,1p2E11.3,0p, 2f8.4,f10.6,1pE11.3)')
     *        buf(i).ldep,  buf(i).code,  buf(i).subcode,
     *        buf(i).charge, buf(i).ridx, buf(i).faiidx,
     *        buf(i).rinmu, buf(i).fai, buf(i).Ek,
     *        buf(i).t, buf(i).wx, buf(i).wy, buf(i).wz,
     *        buf(i).wgt 
#endif
         enddo
      enddo
 100  continue
      end
