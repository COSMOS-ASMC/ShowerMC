!     *********************
!        For the organic and most scintillator, the energy loss by
!      heavy particles (or slow particles?) is not converted to
!      photons  as efficiently as the one by relativistic electrons.
!      This subroutine gives a correction factor for the energy loss,
!      so that you can get effective energy loss by multiplying
!      cf by the true energy loss. 
!         NOTE: this is used for non organic scinti too.
!         called only if quenching terms are given in media data
!        or modifier is specified 
      subroutine epOrgCorrec(modi, media, aPtcl, dedx, cf)
#if defined (Solaris) || (jaxa) || (jaxaflat)
      use epModify,idmodiy=>id
#else 
      use epModify
#endif

      implicit none
#include "Zmedia.h"
#include "Zptcl.h"
#include "Zmass.h"
      integer,intent(in):: modi ! modifier index of the component
                 ! where current ptcl is.
       type(epmedia)::  media  ! input.
       type(ptcl)::  aPtcl   ! input.
      real*8  dedx   ! input.  dE/dx (GeV/(g/cm^2) for the partcle
      real*8  cf     ! output. correction factor. dE_eff = cf x dE_true
      character*1 id
      real(8)::  c1, xx
      real(8)::c2=0. 
      real(8)::cc=0. 
      call epQuenchCoeff(modi, media, aPtcl, dedx, c1, c2, cc,id)
      if( id == "n") then
         cf = 1.
      elseif( id == "T" ) then
         cf = (1.-c2)/(1.+(1.-c2)*c1*dedx) + c2
      elseif( id == "L") then
         xx = c1*dedx + 1.
         cf = xx**( -c2*log(cc*xx))
      elseif( id == "B") then
         cf = 1./(1. + c1 * dedx)
      else
!           should not come but for safety
         write(0,*)' modifyid =', modify(modi)%q%id, 'undef '
         write(0,*) ' ModifyFile index =', modi
         stop 2222
      endif
      end
      subroutine epQuenchCoeff(modi, media, aPtcl, dedx, a,b,c,id)
#if defined (Solaris) || (jaxa) || (jaxaflat)  
      use epModify,idmodiy=>id
#else 
      use epModify
#endif
      implicit none
#include "Zmedia.h"
#include "Zptcl.h"
#include "Zmass.h"
      integer,intent(in):: modi ! modifier index of the component
                 ! where current ptcl is.
       type(epmedia)::  media  ! input.
       type(ptcl)::  aPtcl   ! input. see next
      real*8  dedx   ! input.  dE/dx (GeV/(g/cm^2) for the partcle
               !  these two are not yet used. 
               ! a,b,c may be dependent on these.
      real*8  a,b,c  ! output.  quench coef.
      character*1 id ! output.  one of T,B,L,n for Talre, Birks, LOg
!                               n is for no quenching

      if( modi > 0 ) then
         if( modi > maxModifyNum) then
            write(0,*) 'modifier =', modi, 'for media =',media%name
            write(0,*) ' > ', maxModifyNum, ' i%e. max # in ModifyFile'
            write(0,*) 'To run prg., correct data or give blank '
            write(0,*) ' to ModifyFile in epicsfile'
            stop
         endif
         if( .not. allocated ( modify )) then
            write(0,*) 'Modifier in a component is spcified but '
            write(0,*) '"modify" array has not be allocated'
            write(0,*) 'modify index=', modi
            stop
         endif
         if( IBITS(modify( modi )%kind, bitQuench, 1) > 0  ) then
            if( modify(modi)%q%id == "T") then
               a = modify(modi)%q%a
               b = modify(modi)%q%b
               id = "T"
            elseif( modify(modi)%q%id == "L") then
               a = modify(modi)%q%a
               b = modify(modi)%q%b
               c =  modify(modi)%q%c
               id = "L"
            elseif( modify(modi)%q%id == "B") then
               a = modify(modi)%q%a
               id = "B"
            else
               write(0,*)' modifyid =', modify(modi)%q%id, 'undef '
               write(0,*) ' ModifyFile index =', modi
               stop 2222
            endif
         elseif( media%Birks /= ' ') then  ! not nec. Birks
            call epdefaultQuenchCoef(media, a,b,c,id)
         else
            id = "n"
         endif
      else
         if( media%Birks /= ' ') then
            call epdefaultQuenchCoef(media,  a,b,c,id)
         else
            id = "n"
         endif
      endif
      end
      subroutine epdefaultQuenchCoef(media, a,b,c,id)
      implicit none
#include "Zmedia.h"
       type(epmedia)::  media  ! input.
      real(8),intent(out)::a,b,c  !   quench coef.
      character(1),intent(out)::id !one of T,B,L,n for Talre, Birks, LOg
             ! n is for no quenching
!             use defualt for this media
      a = media%BirksC1
      if( media%Birks == "B" ) then 
!            Birks formula without Chou's 2nd term(c2) or heavy correction
!            term(cc )
!      if( abs( aPtcl.charge ) .gt. 1)  then
!        c1 = c1 * media.BirksCC
!      endif
!      cf = 1./(1. + c1 * dedx + media.BirksC2 * dedx**2)
!        at present, we don't use above two factors
! 
         id = "B"
      elseif(media%Birks == "T" ) then 
         b = media%BirksC2
!               Tarle's formla 
         id = "T"
      elseif( media%Birks == "L" ) then 
!                typcal values c1,c2,cc; 23.53 0.0868 4.611
         b = media%BirksC2
         c = media%BirksCC
         id ="L"
      else
         write(0,*) ' media%Birks =', media%Birks, ' undef'
         stop 1111
      endif
      end
