!     ****************************************************************
!     *
!     * epphotoEp: Photo electric effect prob. per r.l
!     * epphotoEe: samples energy of Photo-electron  and angle
!     *
!     ****************************************************************
!
!
!
!       subroutine epphotoEp(pe, eg, prob, path)
       subroutine epphotoEp(media, eg, prob, path)
       implicit none
#include "Zmedia.h"
#include "Zmass.h"

!       record /photoE/pe  ! input. p.e const
       type(epmedia)::  media   ! from v.8.00
       intent(in):: media

       real*8 eg    ! input. gamma energy in GeV 
       real*8 prob  ! output. prob. per r.l 
       real*8 path  ! output. sampled path in r.l 

       real*8  eout  ! output.  energy of photo electron in GeV
       real*8  cont  ! output.  cosangle of ejected photon relative
                     !  to the incident photon
       real*8 ep, tp, rEe, rEg, a
       integer i
!

       real*8 epi, u, cost



!

       ep=eg/masele

!       if(pe.a + eg .le. masele ) then
!          prob = 1.d-35
!       else
          epi=1./ep
          tp= ( (media%pe%b2*epi + media%pe%b1)*epi + 
     *      media%pe%b0 ) *
     *      media%pe%fa* media%pe%p * epi
!//////////////
!          write(0,*) ' epphotoEp ep=',ep
!          write(0,*) ' eg=',eg, ' tp =',tp
!          write(0,*) ' media.pe.b2=',media.pe.b2
!          write(0,*) ' media.pe.b1=',media.pe.b1
!          write(0,*) ' media.pe.b0=',media.pe.b0
!          write(0,*) ' media.pe.fa=',media.pe.fa
!          write(0,*) ' media.pe.p=',media.pe.p
!          write(0,*) ' media.pe.l=',media.pe.l
!          write(0,*) ' media.pe.cr=',media.pe.cr
!          write(0,*) ' media.Zeff=',media.Zeff
!c////////////////

          if(eg .lt. media%pe%ek ) then
             tp=tp/media%pe%l
          endif
          prob = tp*media%pe%cr
!       endif
       if(media%Zeff .gt. 65.) then
          if(eg .lt.1.e-3) then
!              //////////// this is for Pb
             prob = prob*(1. + (1.e-3/eg)*0.07)
          endif
       endif
       call rndc(u)
       path = -log(u) / prob

       return

!      ************
       entry epphotoEe(media,  eg, eout, cost)
!      ************
!
       if(media%xcom%size .eq. 0 ) then
!            no xcom data
          eout=eg + media%pe%a
          if( eout .le. masele ) then
!              assume eg is transferred to electron
             eout = masele + eg
          endif
       else
          do i = media%pe%noOfShells, 1, -1
             if(eg .gt. media%pe%shellE(i)) then
                eout = eg-media%pe%shellE(i) + masele
!                  if eg=shellE incidentally, eout=masele
!                  and Xray energy will be the same as current eg
!                  then, the loop may happne. to avoid such
!                  we put eg to eout
                if( abs( eout- masele) .lt. 1.d-6) then
                   eout = masele + eg
                endif
                goto 100
             endif
          enddo
!              with this energy and xcom table region, p.e effect
!              cannot take place so eg is very small so we
!              assume Eg is transferred to electron;  
          eout = masele + eg
!c////////////
!          write(0,*) ' ********* eg=',eg
!///////////////
 100      continue
       endif
       rEg=eg/masele
       rEe=(eout-masele) /masele
       if(rEe .le. 0.) then
          cost=1.
       else
          a = ( media%Z2eff/137.0/137.0 + 2.0*rEe + rEg**2 )/
     *    (2.0*rEg*sqrt(2.0*rEe) )
          call ksampPEang(a, cost)
       endif
       end
