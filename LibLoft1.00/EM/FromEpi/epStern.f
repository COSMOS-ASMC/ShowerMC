!     ****************************************************************
!     *                                                              *
!     * epShtern:  to compute coefficient which appear in sternheimers *
!     *          dE/dx expression.
!     *        tested                                                 *
!
! /usage/
!
!      See Sternheimer, p.r. b vol.3(1971)3681; bef v.8.80
!        
!    This is updated in v.8.80 
!     
!
!
!      subroutine epStern(w0in, media)
      subroutine epStern( media)
      implicit none
#include "Zmedia.h"
#include "Zcode.h"
#include "Zptcl.h"
#include "Zmass.h"

!      real*8  w0in  !  input. kinetic energy of knock-on electrons >c w0in
!                    !         is not taken into account as energy locss( GeV).

       type(epmedia):: media  ! input. 
!                             media.I, media.Z,media.A,
!                             media.rho, media.gasF
!                             epExpot must have been  called.
                             ! output.  media.sh
!                     

       type(ptcl):: electron
!
       type(sternh):: sh    !  Sterhnheimer's consts is stored here  tempo.
!                           and copied to media.sh later
!
      real*8  hnp, cb, ep, dedxfull
!

      ep = media%I * 1.d9  ! ep in  eV.
      hnp = media%wp*1.d9  ! plasma energy in eV
!             mdeia.wp in GeV should have been fixed
!             in epGetEffZA
!           next two  the same
!      sh.a = 0.1536* (media.Z/media.A) ! 0.1536 was bef. v8.80
      media%sh%a = 0.153536* media%ZbyAeff
!      sh.b is defined diffrently from v8.80;
!      it comes from a part of the term  ln(2mg^ b^2Tm/I^2)
!      = ln(m^2/I^2) + ln(2*g^2b^2Tm/m) 
!      = 2ln(m/I) + ln(2*g^2b^2Tm/m) 
!      = -2ln(I/m)+ ln(2*g^2b^2Tm/m) 
!        (g= gamma , b=beta, m=massele, Tm=max recoil  kinetic E)
!            
!      the first term is sh.b; second term is taken into account
!      at dE/dx calculation.  (for e-; Tm/m= (g-1)/2;
!      for e+; Tm/m = g-1 )
      media%sh%b = 2*log(masele*1.d9/ep)
!
      sh%c = -2*log(ep/hnp) - 1.
!         Cbar is -sh.c 
      cb =- sh%c
      if(media%gasF .eq. 0) then
!               solid /liquid
         if(ep  .lt. 100.) then
            if(cb .lt. 3.681) then
               sh%x0=0.2
            else
               sh%x0=0.326*cb -1.0
            endif
            sh%x1=2.0
         else
            if(cb .lt. 5.215) then ! bef. v8.80 5.211
               sh%x0=0.2
            else
               sh%x0=0.326*cb-1.5
            endif
            sh%x1=3.0
         endif
!            exception
         if( media%name .eq. "LiqH2" ) then
            sh%k = 5.949
            sh%x0 = 0.425
            sh%x1 = 2.0
         else
            sh%k = 3.
         endif
      else
!             gas 1 atm 0degree ; lot of diff. from versios  bef. 8.80
         if(cb .lt. 12.25) then
            sh%x1=4.
            if(cb .ge. 11.5) then
               sh%x0 = 2.0
            elseif(cb .ge. 11.0) then
               sh%x0 = 1.9
            elseif(cb .ge. 10.5) then
               sh%x0 = 1.8
            elseif(cb .ge. 10.0) then
               sh%x0 = 1.7
            else
               sh%x0 = 1.6
            endif
         else
            sh%x1=5.0
            if(cb .lt. 13.804) then            
               sh%x0 = 2.0
            else
               sh%x0=0.326*cb - 2.5
            endif
         endif
!          exceptions
         if( media%name .eq. "H2") then
            sh%k = 4.754
            sh%x0 = 1.837
            sh%x1 = 3.0
         elseif( media%name .eq. "Hegas") then
            sh%k = 3.297
            sh%x0 = 2.191
            sh%x1 = 3.0
         else
            sh%k = 3.   
         endif
      endif
!      media.sh.xa=cb/4.60517          !  2ln(10)
!      sh.sa = (cb-4.60517X0)/((x1-X0)**k
      sh%sa=(cb-4.60517*sh%x0)/(sh%x1-sh%x0)**sh%k
!!
!       some additional const relating to restriced energy loss

!!         to be set in  epLightNewComp-->epSetEmin.f (v9.154)
!!      media.sh.tcut = w0in        ! in GeV 
!!      media.sh.w0 = w0in * 1000.  ! in MeV unit. 
!!      media.sh.wlg0 = log(media.sh.w0)  !! not used at all


      if( media%format .eq. 2) then
         write(0,*)
     *   "(Sternheimer's) Density effect correction consts. "
         write(0,*)
     *    "Normally diff. between the  input values and ones computed"
         write(0,*) " is rather large, but resultant dE/dx is simlar."
         write(0,*)
     *    " We use input values.  However, check if diff. is large"
!         write(0,*) "If diff < 1%, the values will  not be  printed"
         call epchkSternConst("c",  sh%c, media%sh%c)
         call epchkSternConst("x0", sh%x0, media%sh%x0)
         call epchkSternConst("x1", sh%x1, media%sh%x1)
         call epchkSternConst("a",  sh%sa, media%sh%sa)
         call epchkSternConst("k",  sh%k, media%sh%k)
         call epchkSternConst("delta0",  0.d0, media%sh%delta0)         
      else
         media%sh%c = sh%c
         media%sh%x0 = sh%x0
         media%sh%x1 = sh%x1
         media%sh%sa = sh%sa
         media%sh%k = sh%k
         media%sh%delta0 = 0.
      endif
!         next one is used only when EdepdEdx = f.
!         for compararison with analytical calc.
!         In that case, w0in must be as large as
!         incident particle energy.
      call cmkptc(kelec, regptcl, -1, electron)
      electron%fm%p(4)= 3.162 * masele  
!!!!!!!!!
      media%sh%tcut = 1.0d10
      media%sh%w0 = 1.0d10
!!!!!!!!!!
      call epdedxe(media, electron, media%dEdxatp3m, dedxfull)
      end
!         ***************************
      subroutine  epchkSternConst(name, cmp, inp)
      implicit none
!         check Sternheimer's const computed vs input
!
      character*(*)  name  ! input. variable name
      real*8  cmp  ! input.  computed value in epStern
      real*8  inp  ! input.  corresponding value given by input (fmt=2)
                   !    inp and cmp is compared if the diff 
                   ! output.  if this is -100. cmp is put.

      if( inp .eq. -100.0) then
         inp = cmp
      else
         if(name .ne. "delta0" ) then
!            if(abs(inp-cmp)/inp .gt. 0.01) then
               write(0,*) name, ": input=", inp, " computed=",cmp
!            endif
         else
!            delta0 is always 0 in epStern but expected to be
!            a max of 0.14  for conductors
!
            if( inp .gt. 0.14) then
               write(0,*) "delta0 for ", name, " is > 0.14" 
               write(0,*) "check the value in Media file"
            endif
         endif
      endif               
      end
