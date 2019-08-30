!       treat interaction of  non e,g projetile : pj
!
!     Generated  particles are stored in Pwork(1:Nproduced) which is
!     in modColInfo.  Pwork is type(ptcl)
!     Momentum in Pwork is given in the system where pj is defined.
!     In the case of hadron collisions, particle gneration is
!     done by taking the direcion of projectile as the z-axis.

!     After that, the interface routine rotates the axis so that
!     gnerated particle momentum is seen in the sytem where pj is
!     define. This transormation is done by using crot3mom in the
!     interface routine.   (e.g,  In ghe case of dpmjet,
!     Import/DPM/Interface/cdpmjet.f).
!
!     In the case of decay, such transformation is taken into account
!     by using Lorenz transformation of particles of which momenumt is
!     defeined in the rest system of decaying particle.  And we can
!     include polariazation effect in the momentum of Pwork.
!
      subroutine cinteNEP( pj )
      use modSetIntInf
      implicit none

#include  "Zcode.h"
#include  "Zptcl.h"

!                 
      type(ptcl), intent(inout):: pj ! normmally only "in"  but  in very rare case,
! input partcle is forced to be different one (say collision of Dmeson;
!     its collision cannot be treated yet, so we force to change it to kaon)
      
!

!
      character*80 msg
!

!////////////
!      logical deb
!      common /cccdeb/deb
!/////////////

      
!/////////
!      if(deb) then
!         write(*,*) 'NEP:proc=',  IntInfArray(ProcessNo).process
!      endif
!////////////
!
      if(IntInfArray(ProcessNo)%process .eq. 'knoc') then
         call epNEPKnock( pj ) 
      else
         if(pj%code .eq. kpion) then
            call cintePion(pj)
         elseif(pj%code .eq. kkaon) then
            call cinteKaon(pj)
         elseif(pj%code .eq. knuc) then
            call cinteNuc(pj)
         elseif(pj%code .eq. kmuon) then
            call cinteMuon(pj)
         elseif( pj%code .eq. kgnuc ) then
            call cinteHeavy(pj)
         elseif(pj%code .ge. kalfa .and.
     *           pj%code .le. khvymax  ) then
!cc            elseif(MovedTrack.p.code .ge. kdeut .and.
            call cinteHeavy(pj)
         elseif(pj%code .eq. ktriton) then
            call cinteHeavy(pj)
         elseif(pj%code .eq. kdmes ) then
            call cinteDmes(pj)
         elseif(pj%code .eq. knnb) then
            call cintennb(pj)
         elseif(pj%code .eq. kddb) then
            call cinteddb(pj)
         elseif(pj%code .eq. keta) then
            call cinteEta(pj)
         elseif(pj%code .eq. kgzai) then
            call cinteGzai(pj)
         elseif(pj%code .eq. kbomega) then
            call cinteBomega(pj)
         elseif(pj%code .eq. klambda) then
            call cinteLambda(pj)
         elseif(pj%code .eq. klambdac) then
            call cinteLambdac(pj)
         elseif(pj%code .eq. ksigma) then
            call cinteSigma(pj)
         elseif(pj%code .eq. krho) then
            call cinterho(pj)
         elseif(pj%code .eq. komega) then
            call cinteomega(pj)
         elseif( pj%code .eq. kphi) then
            call cintephi(pj)
         elseif( pj%code  .eq. kds ) then
            call cinteds
         elseif( pj%code .eq. ketap ) then
            call cinteEtap
         elseif( pj%code .eq. kDelta ) then
            call cinteDelta
         elseif( pj%code .eq. ktau ) then
            call cinteTau(pj)
         else
            write(0,*) ' E=', pj%fm%p(4)
            write(msg, *) ' ptcl =', pj%code,
     *        ' interaction=', 
     *         IntInfArray(ProcessNo)%process,
     *        ' should not occure'
            call cerrorMsg(msg, 0)
         endif
      endif
      end     
