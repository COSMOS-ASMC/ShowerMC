      program main
!           test neutron elastic scattring
!          
      implicit none
#include "Zptcl.h"
#include "Zmass.h"
      record /ptcl/ a(2)
      integer ityp, ktyp
      real(8):: ata, atz, epin
      integer:: i, j
      integer:: n
      real(8):: px, py, pz
        real(8):: tetacm, tetalab ! angle in CM and lab
        real(8)::  pt
        real(8):: g, b, s, mp, mt, pin, eincms, pzcms

      ityp=2
      ktyp=2112
      ata=28
      atz =14

      epin = 21.7               ! MeV

      mt = ata*masn*1000. 
      mp = masn*1000.
      pin =sqrt( (epin+mp)**2 - mp**2 )
      s = 2*mt*(epin+mp)+  mt**2 + mp**2
      g = (epin + mp+ mt)/sqrt(s)
      b = sqrt(1.-1./g/g)
      eincms = g*(epin+ mp -b* pin)
      write(0,*) ' g=', g,  ' b=',b, ' Ep cms=', eincms

      call cprePhits
      do i = 1, 100000
         call nelst(ityp,ktyp,ata,atz,epin)
         call cphitsOut(n, a)
         if(n >2 ) then
            write(0,*) ' n>2 strange'
            stop
         endif
         do j = 1, n
c            write(0,*) i, a(j).code, a(j).charge, a(j).fm.p(4)
            if( a(j).code == 6 .and. a(j).charge == 0 ) then
               if(j /= 1 ) then
                  write(0,*) ' scatterd n ? for j=',j
               endif
               pt = sqrt( a(j).fm.p(1)**2 + a(j).fm.p(2)**2)
               tetalab = atan2( pt, a(j).fm.p(3) )*180/3.1415
               pzcms = g*(a(j).fm.p(3) - b*a(j).fm.p(4))
               tetacm =  atan2( pt, pzcms )*180/3.1415
               write(*,'(1p,3g12.3)')
     *         tetalab, tetacm, a(j).fm.p(4)-a(j).mass
            endif
         enddo
      enddo
      end
