#ifdef IBMAIX
        subroutine csmpColInA2(pj, ia,  nc)
        implicit none
#include  "Zptcl.h"
        type(ptcl):: pj  !input projectile particle
        integer ia       ! target mass no.
        integer nc       ! output number of collisions.
        call cerrorMsg('*** IBM AIX must use SucInt=0 ***', 0)
        end
#else

!        ***********************************************************
!        *
!        * csmpColInA2: sample # of collisions inside nucleus
!        *
!        *************** tested 88.08.03***********************k.k**
!                        revised 96.10.16
!           revision:  larger table with fine mesh based on more 
!              accurate computation
!              Air is treated separtely so that p of 10^22 ev can be input.
!
!    Hadron nucleus collision is decompsed into successive
!         collision of incident hadron (p, pi, etc)
!         with nucleon inside the nucleus. This program
!         obtains the # of successive collisions.
!   /usage/ call csmpColInA2(proj, ia, nc)
!       proj: /ptcl/   input. projectile ptcl
!         nc: output.  # of collistions sampled
!   /method/
!       Using wood-saxon density of nucleus, simplified
!         glauber calculation is done by using cwoodsaxon_den etc.
!         its results for 
!         A**(1/3) = 4**(1/3) to 208**(1/3) with step (tatal width)/40
!         and for elementary cross sections log10(15mb) to log10(80mb)
!         step (total width)/20  is tabulated. (cumProb2.h)
!         For air we use separate table for A=14, 16, 40 with
!         log10(15mb)-log10(320mb) / 45
! 
        subroutine csmpColInA2(pj, ia,  nc)
        implicit none

!  #include  "Zcode.h"
#include  "Zptcl.h"
        integer i
        type(ptcl):: pj  !input projectile particle
        integer ia       ! target mass no.
        integer nc       ! output number of collisions.

        real*8 a1, a2, a3, da
        real*8  a, xs
!
        integer  idxa, idxxs

        real*8 u
        real*8  xs1, xs2, dxs
        integer mm, nn, kk
        parameter( mm = 20, nn = 21,  kk= 41)
!                 cumProb2,  xsec     A

        real*4 cumProb2(mm, nn, kk), cumProb14(14, 46),
     *    cumProb16(16,46),   cumProb40(21, 46) 




!        parameter (a1 = 4.0**0.3333333, a2 = 208.**0.333333333333,
        parameter (a1 = 1.5874011, a2 = 5.9249921,
     *  da = (a2-a1)/(kk-1) )

!         parameter ( xs1 =log10(15.), xs2 =log10(80.) )
        parameter ( xs1 = 1.176091259d0, xs2 =  1.903089987d0,
     *  dxs =( xs2 - xs1)/(nn-1))

#include  "cumProb2.h"
#include  "cumProb14.h"
#include  "cumProb16.h"
#include  "cumProb40.h"
!
      if(ia .eq. 14) then
         call csampCollN(pj, cumProb14, 14, 46, nc)
      elseif(ia .eq.  16) then
         call csampCollN(pj, cumProb16, 16, 46, nc)
      elseif(ia .eq. 40) then
         call csampCollN(pj, cumProb40, 21, 46, nc)
      else
!             get cross-section for proton target.        
!        call cxpXsec(pj, xs)
        call cinelx(pj, 1.d0,  1.d0,  xs)
        xs = log10(xs)
!
         a = ia
         a3 = a**0.3333333333
         idxa = (a3- a1)/da + 1
         idxxs = (xs - xs1)/dxs + 1
         if( (a3 - idxa * da - a1) .gt. (idxa*da + da + a1 - a3) ) then
            idxa = idxa  + 1
         endif
         idxa =max(1, min(idxa, kk))

         if( (xs - idxxs * dxs - xs1) .gt.
     *        (idxxs * dxs + dxs + xs1 - xs) ) then
            idxxs = idxxs +1
         endif
         idxxs =max(1, min(idxxs, nn))
         call rndc(u)
         do   i=1,  mm
            if(u .le. cumProb2(i, idxxs, idxa) ) then
               nc=i
               goto 100
            endif
         enddo
 100     continue
      endif
      nc = min(ia, nc)
      end
#endif
      subroutine csampCollN(pj, prob,  np, nx, nc)
      implicit none

#include  "Zptcl.h"

        type(ptcl):: pj  !input projectile particle
        integer np       ! input table size for prob.
        integer nx       ! input. table size for xsection.
        integer nc       ! output number of collisions.
        real*4  prob(np, nx)  ! table for n-coll. prob.

        real*8  xs
!
        integer   idxxs, i

        real*8 u
        real*8  xs1, xs2, dxs
        integer nn
        parameter (nn = 46)  ! this should be nx
!         parameter ( xs1 =log10(15.), xs2 =log10(320.) )
        parameter ( xs1 = 1.176091259d0, xs2 =  2.505149978d0,
     *  dxs =( xs2 - xs1)/(nn-1))
!
        if(nx .ne. nn) then
           write(*,*) ' error input to csampCollN'
           stop 9999
        endif
!             get cross-section for proton target.        
!        call cxpXsec(pj, xs)
        call cinelx(pj, 1.d0, 1.d0,  xs)
        xs = log10(xs)
!
        idxxs = (xs - xs1)/dxs + 1
        if( (xs - idxxs * dxs - xs1) .gt.
     *        (idxxs * dxs + dxs + xs1 - xs) ) then
            idxxs = idxxs +1
         endif
         idxxs =max(1, min(idxxs, nn))

         call rndc(u)
         do   i=1,  np
            if(u .le. prob(i,  idxxs) ) then
               nc=i
               goto 100
            endif
         enddo
 100     continue
         end




