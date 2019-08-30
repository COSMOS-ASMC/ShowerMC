! 
!  gives mean -de/dx  (gev/(g/cm2)) for heavy particle (A,Z>1)
! 
!      Use of SRIM data: 
!           In what follows, Z is charge of the heavy.
!      At config read time, 
!       epSrimChk is called for each media;
! epSrimChk:
!         Make media.srim = -1
!         existence of the directory  Data/Media/'media.name'_srim
!         is examined. (-->see if Z=2 data exist or not;
!                       if data for Z=2 is absent, we regards no data for
!                       srim)
!         If exists, give a sequence # to media.srim
!         and if SRIM data for charge Z exist, 
!            srimno(media.srim).c(Z).size = 0
!         if not
!            srimno(media.srim).c(Z).size=  -1
!         At this time, we don't read data yet until actual request
!         for dE/dx for heavy is made.
!        
! In this subroutine
!            If srim data exits  and its use is not 
!            intentionally avoided (if StoppingPw < 0, or
!            Ek/n > SrimEmax,  srim is
!            not used even if it exits), use Srim data
!            else use old method ; routine for non Srim data

      subroutine epdedxhvy(media, aPtcl, dedt, dedtfull)
      use srimdata
      use modEMcontrol
      implicit none
#include "ZepMaxdef.h"
#include "ZepManager.h"
#include "Zmedia.h"
#include "Zptcl.h"
#include "Zcode.h"
#include "Zmass.h"
!  #include "ZepTrackp.h"


       type(epmedia):: media       ! input. media.sh

       type(ptcl):: aPtcl        ! input. a particle


      real*8 dedt     ! output. restricted <Energy loss rate>  GeV / (g/cm2). 
      real*8 dedtfull ! output.  full dE/ct  
      real*8  kepn, dedtproton, ekt
      real*8  epdedxqeff, temp

      character*80 msg
      logical usesrim 
      integer Z, idx
      real(8):: dedtPfull

       type(ptcl)::  proton !  with E/A

      if(aPtcl%code .ne. kgnuc) then
         write(msg, *) ' ptcl code=',aPtcl%code,
     *     ' is not a heavy particle'
         call cerrorMsg(msg, 0)
      endif
      ekt =(aPtcl%fm%p(4)-aPtcl%mass) ! total K.E
      kepn = ekt/aPtcl%subcode  ! T/A

      proton = aPtcl
      call cmkptc(knuc, -1, 1, proton)
      proton%fm%p(4) = kepn + proton%mass
      call epdedxNone(media, proton, dedtproton, dedtPfull)
!      if(kepn .gt. 0.7) then
!            dedt has Z**2 enhancement already.
!      else
      Z = aPtcl%charge
      if(Z> MAXHEAVYCHG ) then
         usesrim = .false.
      elseif(media%srim <= 0 .or. StoppingPw < 0 ) then
         usesrim = .false.
      elseif(srimno(media%srim)%c(Z)%size == 0 ) then
         call epSrimRead( iowk,  srimno(media%srim), Z)
!         srim data is read for charge Z and make
!         srimno(media.srim)%c(Z)%size=actual size
         idx =srimno(media%srim)%c(Z)%size
!         usesrim= srimno(media.srim)%c(Z)%Ekt(idx) >= ekt
         usesrim = srimno(media%srim)%c(Z)%Ekt(idx) >= ekt 
     *      .and. kepn < SrimEmax   ! kepn > SrimEmax, restricted 
                                    ! loss must be used
      else
         usesrim = srimno(media%srim)%c(Z)%size > 0
         if(usesrim) then
            idx =srimno(media%srim)%c(Z)%size
            usesrim= srimno(media%srim)%c(Z)%Ekt(idx) >= ekt
     *      .and. kepn < SrimEmax
         endif
      endif
      if( usesrim ) then
         call epSrimdEdx(srimno(media%srim), 
     *        aPtcl%fm%p(4)-aPtcl%mass, Z,  dedt)
!              this should be full dedt
         dedtfull = dedt
         dedt = dedt                 
     *    - (dedtPfull - dedtproton) *(Z * epdedxqeff(aPtcl))**2
          !  at high E   ( Z*effZ)**2 ~t Z^2 ;==> restricted dE/dx
          !  at low  E where Srim is important.  
          !     dedtPfull - dedtproton ~ 0 ;==> dedt = dedtfull
                             
      else             
!           various correction for the low energy heavy.
         temp = (Z * epdedxqeff(aPtcl))**2
         dedt = dedtproton * temp
         dedtfull = dedtPfull * temp
      endif

      end
!     *********************************
      real*8 function epdedxqeff2(aPtcl)
!           fractional effective charge 
!           Z*this --> effective charge
      implicit none
#include "Zptcl.h"
       type(ptcl):: aPtcl        ! input. a particle

      real*8  gi, beta, x,  z23, qeff

!               get beta
      gi=aPtcl%mass/aPtcl%fm%p(4)
      if(gi .lt. 1.d-2) then
         beta = 1.d0 - gi**2/2.0
      else
         beta =sqrt(1.d0 - gi**2)
      endif

      z23 = aPtcl%charge**(0.666666)
      x = 121.414*beta/z23 + 0.0378*sin(190.72*beta/z23)
      if(aPtcl%charge == 2) then
         qeff =  (1. -
     *        (1.034- 0.1777*exp(-0.0811*aPtcl%charge))*
     *        exp(-x*1.40) )
      else
         qeff =  (1. -
     *        (1.034- 0.1777*exp(-0.0811*aPtcl%charge))*
     *        exp(-x) )
      endif
      epdedxqeff2 = qeff
      end
!     *********************************
      real*8 function epdedxqeff1(aPtcl)
!          
!         This is from  Pierce and Blann Phys. Rev. 1968 vol.173,No2.
!        pp.390-404.  (fractional effective charge; gamma as they say)
!       ERRORTUM:  (1-exp(-2.5Vp)) should read (1-exp(-5.0Vp))
!          (in abstract and Eq.(6))
!       The paper says
!        Vp must be considered at  E/A<0.3MeV;   Vp must be in sqrt(MeV) unit 
!        I don't know how to use this unit. If we put 1/2 Mp Vp^2 =0.3MeV
!        and we make Mp=1, Vp= sqrt( 0.6MeV) = 0.75sqrt(MeV). 
!        Then (1- exp(-5Vp))=0.976

      implicit none
#include "Zptcl.h"
       type(ptcl):: aPtcl        ! input. a particle

      real*8  gi, beta, x,  z23

      real*8  Vr, Vp, kepn, effZ

!            compute effective charge         
!               get beta
      gi=aPtcl%mass/aPtcl%fm%p(4)
      if(gi .lt. 1.d-2) then
         beta = 1.d0 - gi**2/2.0
      else
         beta =sqrt(1.d0 - gi**2)
      endif

      z23 = aPtcl%charge**(0.666666)
!           Vr = V/(V0*Z23)  
!        V: ion velocity, V0Z23; Tomas-Fermi electron velocity
!        V0 = e^2/hbar (cgs) = e^2 c/hbar*c  = c/137
!        V/C = beta; then Vr = beta*137/z23 
      Vr = beta*137/z23
      if( aPtcl%charge == 2 ) then
!           special treatment for He. 
         effZ =  (1.0 - exp(-1.4*Vr) )
      else         
         effZ  =  (1.0 - exp(-0.95*Vr) )
      endif

      kepn = ( aPtcl%fm%p(4) - aPtcl%mass) /aPtcl%subcode
      if(kepn < 0.5e-3) then
         Vp = sqrt(kepn*1.e3*2)  ! kepn--> MeV
         effZ =  effZ/(1.0 - exp(-5.0*Vp))   
      endif
      epdedxqeff1 =  effZ
!//////////////////////
!      write(*,'(a, 1p,2g14.4, 2i4)')
!     *   'q ', kepn, effZ, aPtcl.charge, aPtcl.subcode
!        effZ is 1 kepn>10MeV/n
!///////////////
      end
!          generic fractional effective charge
      real*8 function epdedxqeff(aPtcl)
      use modEMcontrol
      implicit none
!#include "ZepMaxdef.h"
!#include "ZepManager.h"
!#include "Zmedia.h"
#include "Zptcl.h"
!#include "Zcode.h"
!#include "Zmass.h"
! #include "ZepTrackp.h"

       type(ptcl)::  aPtcl ! input  heavy particle

      real(8):: epdedxqeff1, epdedxqeff2
      if( abs(StoppingPw) == 1 ) then
         epdedxqeff = epdedxqeff1(aPtcl)
      elseif( abs(StoppingPw) == 2 ) then
         epdedxqeff = epdedxqeff2(aPtcl)
      else
         write(0,*) ' StoppingPw=',StoppingPw, ' invalid'
         stop
      endif
      end





