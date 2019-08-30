c          convert Epos result  pdgcode px py pz E     
c         --> code subcode charge pt eta      
c 
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      integer:: evtno, ntp, diffcode
      integer:: i
      real(8):: roots, flag, p(4)
      record /ptcl/ a(50000)
      integer:: code, subc, charge, kf
      real(8):: y, eta, pt
      integer:: npzm, npzp,  Ncht,  Nchpzm,  Nchpzp
      integer,save::eventc = 0
c          1      6 0.450000E+03  2      
      do while (.true.)
         read(*, *, end=1000) evtno, ntp, roots, flag
         diffcode = flag
         npzm = 0
         npzp =0
         Ncht = 0
         Nchpzm = 0
         Nchpzp  =0
         do i = 1, ntp
            read(*, *, end=9000)  kf, a(i).fm.p(:)
            call ckf2cos(kf, code, subc, charge)
            call cmkptc( code, subc, charge, a(i))
            if(charge /= 0) then     
               Ncht = Ncht +1
               if(a(i).fm.p(3) < 0. ) then
                  Nchpzm = Nchpzm + 1
               else
                  Nchpzp = Nchpzp + 1
               endif
            endif
         enddo
         eventc = eventc +1
         write(*,'("h ",i3, 6i6)' )
     *        diffcode, ntp, npzm, npzp, Ncht, Nchpzm, Nchpzp
         do i = 1, ntp
            call cyeta(a(i), y, eta)
            pt =sqrt( a(i).fm.p(1)**2 + a(i).fm.p(2)**2 )
            write(*,'(2i4, 1p, g13.4, 0p, f10.4)')
     *           a(i).code, a(i).charge, pt, eta
         enddo
         write(*,*) 
      enddo
 1000  continue
      write(0,*) ' no.of events =',eventc
      stop
 9000  continue
      write(0,*) ' event next to ',eventc, ' is strange'
      end
