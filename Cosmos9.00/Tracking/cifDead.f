!            cifDead;  see if too low energy 
!
      subroutine cifDead(a, icon)
      use modEMcontrol
        implicit none
#include  "Ztrack.h"
! #include  "Zmagfield.h"        
#include  "Zcode.h"
#include  "Ztrackv.h"
#include  "Ztrackp.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zincidentv.h"

!
!       inform total eneegy to the user when the  a particle
!       kinetic energy < Emin. (chaned)
      type(track)::a   ! input.  a track to be examined
      integer icon !   output.  0 ==> alive  
                   !            1 ==> energy < threshold
                   !            2 ==> path exceeded the limit 
                   !            3 ==> out of upper bound
                   !            4 ==> out of lower bound
                   !            5 ==> angle out of limit 
      real*8  ke, cosfromaxis
      integer jold

      icon = 0
!          for (evaporation) heavy fragment, KE<0 happens
!         due to the fact that p,n mass diff is neglected
!         we don't care.  
      ke = max(a%p%fm%p(4) - a%p%mass, 0.d0)

      if(a%t .gt. PathLimit ) then
         call chookEabsorbD(a, ke, 2)
         icon = 2
      elseif(a%pos%height .ge.  BorderHeightH) then
         call chookEabsorbB(a, 1)
         icon = 3
      elseif(a%pos%height .le.  BorderHeightL) then
         call chookEabsorbB(a, 3)
         icon = 4
      elseif( BackAngLimit .gt. -1.0) then
         if( a%p%fm%p(4) .gt. a%p%mass)  then  ! check only moving one
            call cscalerProd(a%vec%w, DcAtObsXyz, cosfromaxis)
            if(cosfromaxis .lt. BackAngLimit) then
!     if(Eabsorb(1) .ne. 0) then
               if(Eabsorb /= 0) then
                  call chookEabsorbD(a, ke, 4)
               endif
!                discard it 
               icon = 5
            endif
         endif
      endif

      if( icon .ne. 0 )  return    !**************

      if(a%p%code .eq. kphoton .or.  a%p%code .eq. kelec ) then
         if( ke  .lt. KEminCas2) then
            call rndsw(jold, 2)
         endif
         if(a%p%charge .eq. 1) then
            if(a%p%fm%p(4)+a%p%mass .lt. KEminCas) then
!            if(ke .lt. KEminCas) then
               icon = 1
!                     if(btest(Eabsorb(1), BitPositron-1)) then
               if(btest(Eabsorb, BitPositron-1)) then
                  call chookEabsorbD(a, a%p%fm%p(4)+a%p%mass, 0)
!                  call chookEabsorbD(a, ke, 0)
               endif
            endif
         elseif( ke  .lt.  KEminCas ) then 
            icon = 1
            if(a%p%code .eq. kphoton) then
!     if(btest(Eabsorb(1), BitPhoton-1)) then
               if(btest(Eabsorb, BitPhoton-1)) then
                  call chookEabsorbD(a, ke, 0)
               endif
            else
!     if(btest(Eabsorb(1), BitElectron-1)) then
               if(btest(Eabsorb, BitElectron-1)) then
                  call chookEabsorbD(a, ke, 0)
               endif
            endif
         endif
      elseif(a%p%code .eq. knuc) then
!/         if( a.p.subcode .eq. antip) then
         if(a%p%charge .eq. -1 .or. 
     *    ( a%p%charge .eq. 0 .and. a%p%subcode .eq. antip )) then
!     can annihilate
            if( a%p%fm%p(4)+a%p%mass .lt. KEmin ) then
!            if( ke .lt. KEmin ) then
               icon = 1
!     if(btest(Eabsorb(1), BitAntiNuc) ) then
               if(btest(Eabsorb, BitAntiNuc-1) ) then
                  call chookEabsorbD(a, a%p%fm%p(4)+a%p%mass, 0)
!                  call chookEabsorbD(a, ke, 0)
               endif
            endif
            if( a%p%fm%p(4) .lt. KEmin2) then
               call rndsw(jold, 2)
            endif
         else
            if( ke .lt.  KEminObs(knuc) ) then
               icon = 1
               if(a%p%charge .eq. 1 ) then
!     if(btest(Eabsorb(1), BitProton) ) then
                  if(btest(Eabsorb, BitProton-1) ) then
                     call chookEabsorbD(a, ke, 0)
                  endif
!     elseif (btest(Eabsorb(1), BitNeutron) ) then
               elseif (btest(Eabsorb, BitNeutron-1) ) then
                  call chookEabsorbD(a, ke, 0)
               endif
            endif
            if(ke .lt. KEmin2) then
               call rndsw(jold, 2)   
            endif
         endif
      elseif(a%p%code .eq. kpion  .or.
     *        a%p%code .eq. kkaon  .or.
     *        a%p%code .eq. kmuon) then
!               can decay;  max energy of  last electron is Et incident
!            allmost all decay before coming here. don't worry 
!              2013/Feb.26
!         if(a.p.fm.p(4) .lt. KEmin) then
         if( a%p%fm%p(4) .lt. min(KEmin,KEminObs(a%p%code) )) then
!         if(ke .lt. KEmin) then
            icon = 1
            if(a%p%fm%p(4) .lt. a%p%mass) then
!             very rare but happened (~2 times in 100 TeV p)
               a%p%fm%p(4) = a%p%mass
            endif
!     if(btest(Eabsorb(1), BitDecay-1)) then
            if(btest(Eabsorb, BitDecay-1)) then
               call chookEabsorbD(a, ke,  0)
!               call chookEabsorbD(a, a.p.fm.p(4),  0)
            endif
         endif
         if(a%p%fm%p(4) .lt. KEmin2) then
            call rndsw(jold, 2)
         endif
      elseif( a%p%code .eq. kgnuc ) then
!              heavey
         if( ke .lt.  KEminObs(knuc)) then
            icon = 1
!     if(btest(Eabsorb(1), BitProton) ) then
            if( btest(Eabsorb, BitProton-1) ) then
               call chookEabsorbD(a, ke, 0)
            endif
            if(ke .lt. KEmin2) then
               call rndsw(jold, 2)   
            endif
         endif
      else         
!                other particles;  basically they should decay
!                so we use here total energy ???
!         if(a.p.fm.p(4) .lt.  KEmin) then
         if(ke .lt.  KEminObs(8)) then
            icon = 1
!     if(btest(Eabsorb(1), BitOther-1)) then
            if(btest(Eabsorb, BitOther-1))then
!               call chookEabsorbD(a, a.p.fm.p(4), 0 )
               call chookEabsorbD(a, ke, 0 )
            endif
         endif
         if(ke .lt. KEmin2) then
            call rndsw(jold, 2)
         endif
      endif

      if( a%where .gt. EndLevel2) then
         call rndsw(jold, 2)
      endif
      end
