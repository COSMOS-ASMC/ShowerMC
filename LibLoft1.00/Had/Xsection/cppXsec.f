!      real*8  e, xs
!      e = 0.1
!      do while (.true.)
!         call cpbarpXsec(e, xs)
!         write(*,*) sqrt( e * (e+2*0.938)), xs
!         e = e * 10.**0.1
!         if(e .gt. 5000.) goto 10
!      enddo
! 10   continue
!      end
!
!   inelastic cross-section at enetire  energy.
!   IncreaseXsec =2 should not be used.
!        cppXsec   p   p      n p
!        cpbarpXsec  p_b p      n_b  n    n_b p   p_b n
!        cpiMinuspXsec  pi- p
!        cpiPluspXsec  pi+ p
!        ckMinuspXsec  k-  p      k0_b n
!        ckPluspXsec k+  p k+ n  k- n   k0_b  p   k0 n
!        real*8 e, xs
!        e=1000.
!        call cppXsec(e, xs)
!        write(*,*) xs
!        end
!       ***********************************************************
!       *
!       * p_p inelastic cross-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!    e:real*8.  input. kinetic energy of partilce in GeV.
!   xs:real*8   output. Inelastic cross-section in mb
!
        subroutine cppXsec(e, xs)
        implicit none
#include  "Zmass.h"
        real*8 e, xs
        real*8 p
!        p = sqrt( (e+masp)**2 - masp**2 )
        p = sqrt( e*(e+masp+masp) )
        call cppInelaXs(p, xs)
      end
!       ***********************************************************
!       *
!       * p_b p inelastic cross-section
!       *
!       ******************                 *****************k.k ***
!
!   e: input.  e kinetic energy in GeV
!   xs:output. inelastic cross-section in mb
!
        subroutine cpbarpXsec(e, xs)
        implicit none
#include "Zmass.h"
        real*8 e, xs
        real*8 p
        p = sqrt( e*(e+masp+masp) )
        call cpbarpInelaXs(p, xs)
        end
!       ***********************************************************
!       *
!       * pi- p inelastic cross-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  pi- k.e     in GeV
!   xs:output. inelastic cross-section in mb
!
        subroutine cpiMinuspXsec(e, xs)
        implicit none
#include "Zmass.h"
        real*8 e, xs
        real*8 p
        p = sqrt( e * ( e + maspic + maspic) ) 
        call cpimpInelaXs(p, xs)
      end
!       ***********************************************************
!       *
!       * pi+ p inelastic cross-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  pi+ k.e     in GeV
!   xs:output. inelastic cross-section in mb
!
        subroutine cpiPluspXsec(e, xs)
        implicit none
#include  "Zmass.h"
        real*8 e, xs
        real*8 p
        p = sqrt( e * ( e + maspic + maspic) ) 
        call cpippInelaXs(p, xs)
      end
!       ***********************************************************
!       *
!       * k-  p inelastic cross-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  k- k.e   GeV
!   xs:output. inelastic cross-section in mb
!
        subroutine ckMinuspXsec(e, xs)
        implicit none
#include  "Zmass.h"
        real*8 e, xs
        real*8 p
        p=sqrt( e*(e+maskc+maskc) )
        call ckmpInelaXs(p, xs)
      end
!       ***********************************************************
!       *
!       * k+  p inelastic cross-section
!       *
!       ****************** tested 88.08.11 *****************k.k ***
!
!   e: input.  k+ k.e    GeV
!   xs:output. inelastic cross-section in mb
!
        subroutine ckPluspXsec(e, xs)
        implicit none
!----        include 'Zxsectionp.h'
#include  "Zmass.h"
        real*8 e, xs
        real*8 p
        p=sqrt( e*(e+maskc+maskc) )
        call ckppInelaXs(p, xs)
      end
