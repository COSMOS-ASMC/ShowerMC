	subroutine epUrban(urb, avededx, dx, aPtcl, de)
	implicit none
#include "Zptcl.h"
#include "Zurban.h"
!    By using Urban model, energy loss is calculated with 
!    full fluctuation; this is consistent with Landau/Vavilov/Gaussina
!    theory if the all delta ray is taken into account.  
!    ***************
!   This must be called after epsetUrban has been called to fix the 
!   medium consts. 
!
!   -----------------------------------------------------------
!   The cut-off energy of delta electron must be
!   lower than the maximum transferable energy to electrons.
!   -----------------------------------------------------------
!
!    ***************
       type(urban):: urb ! input. Urban's const of the medium
	real*8 avededx  ! input. average dE/dx in GeV/cm.

	real*8 dx       ! input. pass length in cm for which energy loss
                        !        is to be sampled.
       type(ptcl):: aPtcl ! input. mass and energy are used.
	real*8 de       ! output. sampled energy loss for pass dx in GeV.
              
	real*8 gamma, beta2, gbeta2
        integer i
        real*8 mgb2
	include "ZdEdxSpec.h"
!
!///////////////
!	write(0,*) ' entering urb: dx =', dx
!	write(0,*) ' code, sub, chg KE=', 
!     *   aPtcl.code,aPtcl.subcode,aPtcl.charge, aPtcl.fm.p(4)-aPtcl.mass
!	write(0,*) ' avededx=',avededx
!///////////////

	if(aPtcl%fm%p(4) - aPtcl%mass .le. 1.d-6 ) then
	   de = avededx*dx
        elseif( dx .le. 0.) then
	   de = 0.
	else
	   if( Tupper  .gt.  urb%EmaxN ) then
	      urb%Emax = urb%EmaxN
	      urb%EmaxI = urb%EmaxIN 
	      urb%IEmaxI = urb%IEmaxIN 
	      urb%EmaxByEmaxI = urb%EmaxByEmaxIN
	      urb%PreSigma3 =  urb%PreSigma3N
	   else
	      urb%Emax = Tupper
	      urb%EmaxI = urb%Emax + urb%AveIoniz
	      urb%IEmaxI = urb%AveIoniz * urb%EmaxI
	      urb%EmaxByEmaxI = urb%Emax/urb%EmaxI
	      urb%PreSigma3 = urb%EmaxByEmaxI /
     *         (urb%AveIoniz * log(urb%EmaxI/urb%AveIoniz) )
     *          *  urb%ExToIoniz
	   endif

	   gamma =aPtcl%fm%p(4)/aPtcl%mass
	   beta2 = 1.0d0 - 1./gamma**2
	   gbeta2 = gamma**2 * beta2
	   mgb2 = aPtcl%mass*gbeta2 * 2

	   do i = 1, 2
	      urb%MacroSigma(i) =
     *        avededx* urb%OscStr(i)/urb%ExcitErg(i) * 
     *        (  log(mgb2/urb%ExcitErg(i)) - beta2 ) /
     *        (  log(mgb2/urb%AveIoniz   ) - beta2 ) *
     *        (1.-urb%ExToIoniz)
	   enddo
	   urb%MacroSigma(3) = avededx * urb%PreSigma3
	   do i = 1, 3
	      urb%AveN(i) = urb%MacroSigma(i) * dx
	      if(urb%AveN(i) .le. 0.) then
		 urb%ElossTotal = avededx*dx
		 goto 100
	      endif
	   enddo	
!            very thin cases must be treated specially
!            by epUrbanThin below, 
	   if(urb%AveN(1)+ urb%AveN(2) + urb%AveN(3) .lt. 4.6) then
!               ! no energy loss probability is > 0.01. treat specially
	      call epUrbanThin(urb, avededx*dx)
	   else
	      call epUrbanNorm(urb)
	   endif
 100	      continue
	   de = urb%ElossTotal
	endif
!////////////
!	write(0,*) ' exiting urban'
!////////////
        end 
!       ------------------------------------------------
        subroutine epUrbanThin(urb, avede)
	implicit none
#include "Zurban.h"
       type(urban):: urb
	real*8 avede  ! input. <dE/dX> *dX  GeV
!          very thin path case.
!          interaction is assumed to be with outer electrons only.
!          Their energy level is 10 eV.
        integer i
	real*8 u
!             original prescription
!	urb.AveNThin = avede/urb.EOuter
!
	urb%AveNThin = avede* (1./urb%EOuter - 1./urb%Emax)
     *    /log(urb%Emax/urb%EOuter)
	call kpoisn(urb%AveNThin, urb%NThin)
	urb%ElossTotal = 0.
	do i = 1, urb%NThin
	   call rndc(u)	
!            original     
!	   urb.ElossTotal = urb.ElossTotal + urb.EOuter/
!     *         (1. - u*urb.Emax/(urb.Emax + urb.EOuter))

	   urb%ElossTotal = urb%ElossTotal + urb%EOuter/
     *         (1. - u*(1.-urb%EOuter/urb%Emax))
        enddo
!           another method; but we need to know peak position
!          we use approximate Lanau distriubution as given
!           by ksplandau
!
!	real*8 b, c
!	c = avede/...
!	b = c*..
!	call ksplandau(b, c, urb.ElossTotal)

        end
!       --------------------------
        subroutine epUrbanNorm(urb)
	implicit none
#include "Zurban.h"
       type(urban):: urb

	integer i
	real*8 u, logAlfa, temp


!          sample actual number of excitation and ionization.
        do i = 1, 3
   	    call kpoisn(urb%AveN(i), urb%Ncoll(i))
        enddo
!         energy loss due to exciteation.
        urb%ElossByExcit = urb%Ncoll(1) * urb%ExcitErg(1)  +  
     *                urb%Ncoll(2) * urb%ExcitErg(2)

!          now loss by ionization.

        if(urb%Ncoll(3) .lt. urb%LargeN3) then
	   urb%ElossByIoniz = 0.
	   do i = 1, urb%Ncoll(3)
	      call rndc(u)
	      urb%ElossByIoniz  = urb%ElossByIoniz +
     *          urb%AveIoniz/(1.0 - u *urb%EmaxByEmaxI)
	   enddo
        else
            urb%Alfa = (urb%Ncoll(3) + urb%Smallc2)*urb%EmaxI/
     *      (urb%Smallc2*urb%EmaxI+ urb%Ncoll(3)*urb%AveIoniz)

            logAlfa = log(urb%Alfa)

            urb%PAtype = urb%EmaxI*(urb%Alfa-1.)/urb%Emax/urb%Alfa

            urb%AveElossAtype = 
     *       urb%AveIoniz *urb%Alfa *logAlfa/(urb%Alfa -1.0)

            urb%SigmaElossA2 = urb%AveIoniz2 * urb%Alfa *
     *        (1. - urb%Alfa* (logAlfa/(urb%Alfa-1.0))**2 )

!             average number of A type col.
            urb%AveAtype = urb%Ncoll(3) * urb%PAtype
            urb%SigmaAtype =
     *      sqrt(max(urb%AveAtype * (1.0d0 - urb%PAtype), 0.d0) )

            call kgauss(urb%AveAtype, urb%SigmaAtype, temp)
	    urb%NAtype =max( temp + 0.5d0, 0.d0)
!             average energy loss by NAtype coll.
	    urb%AveElossByNAtype = urb%NAtype * urb%AveElossAtype
            urb%SigmaByNAtype = sqrt(urb%NAtype * urb%SigmaElossA2)
!              sample from Gaussian
            call kgauss(urb%AveElossByNAtype, urb%SigmaByNAtype, 
     *         urb%ElossByAtype)
	
            urb%ElossByBtype = 0.
	    urb%NBtype = urb%Ncoll(3) - urb%NAtype 
	    do i = 1, urb%NBtype
	       call rndc(u)
	       urb%ElossByBtype = urb%ElossByBtype + 
     *         urb%Alfa*urb%AveIoniz / 
     *         ( 1.0 - u *
     *         (urb%EmaxI - urb%Alfa*urb%AveIoniz)/urb%EmaxI)
            enddo
            urb%ElossByIoniz = urb%ElossByAtype + urb%ElossByBtype
        endif
 
	urb%ElossTotal = urb%ElossByExcit + urb%ElossByIoniz
!       /////////////////
! 	write(12, *)
!     * sngl(urb.ElossByIoniz*1.d3/dx), 
!     * sngl(urb.ElossByAtype*1.d3/dx), 
!     * sngl(urb.ElossByBtype*1.d3/dx),
!     * sngl(urb.ElossByExcit*1.d3/dx),
!       /////////////
        end
        subroutine epsetUrban(media, urb)
!
!          set Urban's consts. must be called after epStern has been called
!
        implicit none
#include "Zmedia.h"
       type(epmedia):: media  ! input.  media
       type(urban):: urb !  outupt. Urban's consts.    = media.urb
!
!       EmaxN, EmaxIN, IEmaxIN, EmaxByEmaxIN, PreSigma3N; these are
!         to replace Emax etc (those without N)
!
        urb%LargeN3 = 16
!	urb.ExToIoniz = 0.4
	urb%ExToIoniz = 0.5
        urb%EOuter = 10.d-9   ! 10 eV for outer electron ionization


!        urb.Zeff = media.Zeff
        urb%Zeff = media%Z
	urb%AveIoniz = media%I
	urb%AveIoniz2 = urb%AveIoniz**2
!       here  epResetUrban equivalent one is placed
!       it was ok when only 1 fixed Tcut is used
!       for each media. but it's danger so next one
!       is moved to epResetUrban and called
!       when a new component appears
!!!!!!!!!!!!  v9.154
!	urb.EmaxN = media.sh.tcut  ! in GeV unit
!	urb.EmaxIN = urb.EmaxN + urb.AveIoniz
!	urb.IEmaxIN = urb.AveIoniz * urb.EmaxIN
!	urb.EmaxByEmaxIN = urb.EmaxN/urb.EmaxIN
!	urb.PreSigma3N = urb.EmaxByEmaxIN /
!     *  (urb.AveIoniz * log(urb.EmaxIN/urb.AveIoniz) ) *
!     *   urb.ExToIoniz
!!!!!!!!!!


	urb%Smallc2 = urb%LargeN3
	if(urb%Zeff .le. 2.) then
	    urb%OscStr(2) = 0.
	else	
	    urb%OscStr(2) = 2./urb%Zeff
        endif
	urb%OscStr(1) = 1. - urb%OscStr(2)
!        urb.ExcitErg(2) = 10.d-9 * media.Z2eff  ! approx K-shell energy in GeV
        urb%ExcitErg(2) = 10.d-9 * media%Z2  !  which is better ?
        urb%ExcitErg(1) = 
     *  (urb%AveIoniz/urb%ExcitErg(2)**urb%OscStr(2))
     *          **(1./urb%OscStr(1))
        urb%ExcitErg(1) = 
     *  (urb%AveIoniz*1.d9/(urb%ExcitErg(2)*1.d9)**urb%OscStr(2))
     *          **(1./urb%OscStr(1))*1.d-9

        end

	subroutine epResetUrban(media, urb)  ! v9.154
!         Urban  consts but changing when recoil
!         min. energy changes. must be called
!         from epLightNewComp->epSetEmin->epSetTcut
!        implicit none
#include "Zmedia.h"

       type(epmedia):: media  ! input.  media
       type(urban):: urb !  outupt. Urban's consts.    = media.urb

	real(8),save:: tcutsave = 1.d10
	character(8),save:: namesave="xxxx" 

	if( tcutsave == media%sh%tcut .and.
     * 	    namesave == media%name ) then
	   return  ! **************
	else
	   namesave = media%name
	   tcutsave = media%sh%tcut
	endif
	urb%EmaxN = media%sh%tcut  ! in GeV unit
	urb%EmaxIN = urb%EmaxN + urb%AveIoniz
	urb%IEmaxIN = urb%AveIoniz * urb%EmaxIN
	urb%EmaxByEmaxIN = urb%EmaxN/urb%EmaxIN
	urb%PreSigma3N = urb%EmaxByEmaxIN /
     *  (urb%AveIoniz * log(urb%EmaxIN/urb%AveIoniz) ) *
     *   urb%ExToIoniz
	end subroutine epResetUrban
      subroutine epWriteUrbCnst(urb)
      implicit none
#include "Zurban.h"
      type(urban),intent(in):: urb
      write(0,*) ' urb const '
      write(0,*) ' OscStr(2)=', urb%OscStr(:)
      write(0,*) ' ExcitErg(2)=', urb%ExcitErg(:)
      write(0,*) ' MacroSigma(3)=',urb%MacroSigma(:)
      write(0,*) ' AveN(3) =', urb%AveN(:)
      write(0,*) ' ExToIoniz=', urb%ExToIoniz
      write(0,*) ' AveIoniz=', urb%AveIoniz
      write(0,*) ' Zeff, smallc2=',   urb%Zeff, urb%Smallc2
      write(0,*) ' Emax, EmaxI=', urb%Emax, urb%EmaxI
      write(0,*) ' IEmaxI=', urb%IEmaxI
      write(0,*) ' EmaxByEmaxI, PreSigma3=',
     *     urb%EmaxByEmaxI, urb%PreSigma3
      write(0,*) ' EmaxN, EmaxIN, IEmaxIN=',
     *     urb%EmaxN, urb%EmaxIN, urb%IEmaxIN
      write(0,*) ' EmaxByEmaxIN, PreSigma3N=',
     *     urb%EmaxByEmaxIN, urb%PreSigma3N
      write(0,*) ' ElossByExcit, ElossByIoniz',
     *     urb%ElossByExcit, urb%ElossByIoniz
      write(0,*) '  Alfa, PAtype, AveIoniz2=',
     *     urb%Alfa, urb%PAtype, urb%AveIoniz2
      write(0,*) ' SigmaElossA2, AveAtype=',
     *     urb%SigmaElossA2, urb%AveAtype
      write(0,*) ' AveElossAtype, SigmaAtype=',
     *     urb%AveElossAtype, urb%SigmaAtype
      write(0,*)
     *     ' SigmaByNAtype,AveElossByNAtype, ElossByAtype=',
     *     urb%SigmaByNAtype, urb%AveElossByNAtype, urb%ElossByAtype
      write(0,*) ' ElossByBtype, ElossTotal=',
     *     urb%ElossByBtype, urb%ElossTotal
      write(0,*) ' AveNThin, EOuter=',
     *    urb%AveNThin, urb%EOuter
      end


