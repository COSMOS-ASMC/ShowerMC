c          convert Epos result  pdgcode px py pz E     
c         -->  compatible to Gencol output
c              code subcode charge KE  wx wy wz  user
c           **** without header ***
c     also see ///// below ; to make output short
c  Usage:  
c     make -f convEposRes2.mk
c     ./convEpos2PCLinuxIFC <  eposresult.data
c
      implicit none
#include "Zptcl.h"
#include "Zcode.h"
      integer:: evtno, ntp, diffcode
      integer:: i
      real(8):: roots, flag, p(4)
      record /ptcl/ a(50000)
      integer:: code, subc, charge, kf
      real(8):: y, eta, pt, pabs
      integer:: npzm, npzp,  Ncht,  Nchpzm,  Nchpzp
      integer,save::eventc = 0
      character*80 line
      character*20 field(4)
      integer acnf
c          1      6 0.450000E+03  2      
      do while (.true.)
         read(*,'(a)',end=1000)  line
         eventc = eventc +1
         call kgetField(line, field, 4, acnf)

         read(field(2),*) ntp
         read(field(3),*) roots
         read(field(4),*) flag
c         read(*, *, end=1000) evtno, ntp, roots, flag
c/////////////
c         if(eventc > 1400000) goto 1000
c////////////
         evtno = eventc
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

c         write(*,'("h ",i3, 6i6)' )
c     *        diffcode, ntp, npzm, npzp, Ncht, Nchpzm, Nchpzp
         do i = 1, ntp
c            call cyeta(a(i), y, eta)
            pt =sqrt( a(i).fm.p(1)**2 + a(i).fm.p(2)**2 )
            pabs = sqrt(pt**2 + a(i).fm.p(3)**2)
c/////////////////
            if(  a(i).fm.p(3)/pabs < 0. )  then
            elseif( a(i).fm.p(4) < 10. ) then
            elseif( a(i).code == 4 .and. a(i).charge /= 0 ) then
            else
c//////////////
               write(*,'(3i4, 1p, g13.4, 3g17.8, i4)')
     *           a(i).code,  a(i).subcode, a(i).charge,
     *           a(i).fm.p(4)- a(i).mass, 
     *           a(i).fm.p(1)/pabs, a(i).fm.p(2)/pabs,
     *           a(i).fm.p(3)/pabs, i
c/////////////
               endif
c//////////////
         enddo
         write(*,*) 
      enddo
 1000  continue
      write(0,*) ' no.of events =',evtno
      stop
 9000  continue
      write(0,*) ' event next to ',eventc, ' is strange'
      end
