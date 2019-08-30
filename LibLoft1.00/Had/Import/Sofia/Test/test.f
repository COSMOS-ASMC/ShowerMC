      program testsofia
      implicit none
#include "Zcode.h"
#include "Zmass.h"
#include "Zptcl.h"
      record /ptcl/work(1000), pj, p0, px
      real(8):: eps
      integer  A, Z, ntp, Imode
      integer nev, nevent,j

      call cmkptc(1,0,0, pj)  ! photon
      
      write(0,*)
     *  " enter Eg(GeV), target A, Z, # of events"
      read(*,*) eps, A, Z, nevent
      write(0,*) ' Eg=',eps, 'A Z=',A, Z, ' # of events=',
     *     nevent
!                    P(1:2) need not be 0.  
      pj.fm.p(1) = 0.
      pj.fm.p(2) = 0.
      pj.fm.p(3) = eps
      pj.fm.p(4) =  eps
      write(*,'(a)')
     *   "# ev#  code subcode chg  Pz/P0 Px Py Pz E mode N"
      do  nev=1, nevent
         call csofia(A, Z, pj, work, ntp)
         call csofiamode(Imode)
!            p0 seen from pj
         call cirot3vec(1, pj.fm.p(1), pj.fm.p(1), p0)

         do j = 1, ntp 
!                px seen from pj
            call cirot3vec(j, pj.fm.p(1), work(j).fm.p(1), 
     *                        px.fm.p(1))

            write(*,'(4i4,1p, 5E12.3,i4,i5)' )
     *         nev, work(j).code, work(j).subcode, work(j).charge,
     *         px.fm.p(3)/p0.fm.p(3),   ! x seen from pj direction
     *         work(j).fm.p(1),
     *         work(j).fm.p(2),
     *         work(j).fm.p(3),
     *         work(j).fm.p(4),  Imode, ntp
         enddo
      enddo
       write(*,'(a)')
     *   "# ev#  code subcode chg  Pz/P0 Px Py Pz E mode N"
       END
