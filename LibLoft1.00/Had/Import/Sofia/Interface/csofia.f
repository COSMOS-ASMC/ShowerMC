!       interface to sophia (photo hadron production routin)
      subroutine csofia(massN, atomicN, pj,  a, ntp)
      use modsofia
      implicit none
#include "Zcode.h"
#include "Zptcl.h"
#include "Zmass.h"
      integer,intent(in):: massN ! A of the target
      integer,intent(in):: atomicN ! Z of the target
      type (ptcl):: pj          ! input.  photon
      type (ptcl):: a(*)        ! output.   produced particles
                           ! they are at the same coord. as pj
      integer,intent(out):: ntp ! # of procuded partilces



!
! E0 = energy of incident proton (in lab frame) [in GeV]
! L0 = code number of the incident nucleon (L0=13: proton, L0=14: neutron)
!     common with   sofia (S_PLIST --> SS_PLIST; why?  see Readme)
      real(8)::P
      integer:: LLIST, NP, Ideb
      COMMON /SS_PLIST/ P(2000,5), LLIST(2000), NP, Ideb


      integer:: ic 
      logical,save:: first=.true.
      real(8),parameter:: theta = 0.
      real(8):: eps  ! photon energy
      real(8):: E0   ! nucleon energy (mass)
      integer:: L0   ! nucleon code
      integer,save::imode
      integer,intent(out):: mode
      integer:: i
!///////////
!      real(8)::sumTE
!//////////

      if( first ) then
         call setDecayCond
         first = .false.
      endif
!          check incident
      if( pj%code /= kphoton ) then
         write(0, *) ' incident code =', pj%code, ' invalid'
         write(0, *) ' for csofia '
         call cerrorMsg(' strage',0)
      endif
!           fix target (p or n)
      call cfxTgtChg(massN, atomicN, ic)
!///////////////
!      write(0,*) 'csofia: pj.E=', pj.fm.p(4)
!      write(0,*) ' tgt ic =',ic
!/////////////////
      

      if(ic == 1) then
         L0 = 13                ! p       
         E0 = masp   ! target is at rest
      else
         L0 = 14   ! n
         E0 = masn
      endif


      call initial(L0)

      eps = pj%fm%p(4)

!***********************************************************
!*** CALL PHOTOPION EVENT GENERATOR 
!***********************************************************
!//////
!      write(0,*) 'calling sofiaEvent=', L0,E0, eps, theta
!///////
      call sofiaEvent(L0, E0, eps, theta, imode)
!///////////
!      write(0,*) ' imode =', imode
!///////////////
!//////////////////
!      sumTE = 0.
!///////////////
      ntp = 0
      do i = 1, NP
         call csibyllcode2cos(LLIST(i), a(i))
         if(a(i)%code /= krare) then
            ntp = ntp + 1
            a(ntp)%fm%p(1) = -P(i,1) ! photon direction is -z; (theta=0)
            a(ntp)%fm%p(2) = -P(i,2) ! to make it +, invert sign
            a(ntp)%fm%p(3) = -P(i,3)
            a(ntp)%fm%p(4) =  P(i,4)
         endif
!/////////////
!         sumTE = sumTE + P(i,4)
!/////////////
       enddo
!c//////////
!       if( abs( sumTE / (eps+E0) - 1.) > 0.005 ) then
!          write(*,'(a,1p,2g13.4, i6)'  )
!     *      'errgp ', eps,  (sumTE / (eps+E0) - 1.), NP 
!       endif
!c/////////
       call crot3mom( pj, a, ntp)

       return
!                 not recommended way but for the moment use entry
       entry csofiaMode(mode)  ! return imode
       mode = imode
       return
       end

