!      module BPKochMotz
      module BPPS   ! partial screening by Messel/Koch/Motz
      implicit none
!       for given media, these quantities are
!       calculated for each element in the media.
!       The calculation is perfomed by Zpart
!       and is used to get brems and pair cross-sections
!       These are independent of energy (and so x=Eg/Ee)
      private
      real(8),save:: Z=0.
      real(8),save::  Z13, Z23, ame, apme, delta, cf,
     *        lnZ43,  lnZ83, Z2,
     *        al183z, ccz, bcoef, x0g, fz
!       x0g:  r.l for virtual matte:r A=1 and Z.
!          next is exception.  dependent on Eg or Ee
      real(8):: Ebyme    ! Ee/Me or  Eg/Me

      public epBrem, epPair, epBpZpart
      contains

      function epBrem(Zin, Eeme, x) result(ans)
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmass.h"

      real(8),intent(in):: Zin ! charge
      real(8),intent(in):: Eeme  ! Ee/me
      real(8),intent(in):: x     ! Eg/Ee
      real(8)::ans
!         
!         first get
!         bremsung function per radiation length.  approximation j in
!         koch & motz notation + correction at low energy
!     Then convert it into mb/target.
!             
!          Since  sigma = alpha r0**2 f(z) g(x) 
!        and 1/X0 = alpha r0**2 f(z)  N0/A (cm2/g)
!     (N0; Abogadro no. A: mass number)
!        and epBrem first gives g(x).  So, sigma can be
!    obtained by  sigma = A/(N0X0) g(x).
!    Since the radiation length has linear A dependence
!       sigma =  g(x)/ N0 / (X0/A); X0/A has no A dependence.
!    x0g  has been given by using A= 1 so that
!    we may devide g(x) by N0*X0 to get sigma in cm^2. 
!    To convert it to mb, 10^{27} is to be multiplied.
!    Hence  10^27/N0 = 1.66e3,  g(x)*1.66e3/x0g is in mb
!
!
!
      real(8)::  v
!      real(9):: bigf, bcorec   ! internal so don't declare here
      real(8),parameter::cmTomb = 1.d27/N0
      real(8):: E

      call epBPZpart(Zin)
      Ebyme = Eeme
      v = x
      E = Eeme * masele
!/////////// without these lines, epBrem becoms 0   Why ???  (in old days)
      if(  bigf(0,v) <=  0. .or.  bcorec(E) <=   0. ) then
         write(0,*)' bigf,bcorec= ',
     *   bigf(0,v),bcorec( E), ' at v,E=', v, E
         stop
      endif
!//////////////   
      ans = max( bigf(0,v)/al183z *bcorec( E )/v, 0.d0)
!          convert it to mb
      ans = ans * cmTomb /x0g
      end function epBrem

!     ************
      function epPair(Zin, Egme, x) result(ans)
!     ************
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"
#include "Zmass.h"

      real(8),intent(in):: Zin  ! target charge
      real(8),intent(in):: Egme  ! Eg/me
      real(8),intent(in):: x     ! Ee/Eg
      real(8):: ans

      real(8),parameter::cmTomb = 1.d27/N0
!
!      real(8), bigg, pcorec  ! interanl so don't declare here
      real(8):: E

      Ebyme = Egme
      E =  Egme * masele        ! Eg not electron E
!        pair creation function  v=e/eg
      ans = max( bigg(x)/al183z* pcorec(E), 0.d0)
!     convert it into mb
      ans = ans * cmTomb /x0g
      end function epPair
!     ****************************************************************
!     *                                                              *
!     *  bigf:  auxliary function for brems cross-section            *
!     *  bigg:  //                    pair                           *
!     *                                                              *
!     ****************************************************************
!
!
!
      function bigf( n, v) result(ans)
      implicit none
      integer,intent(in):: n      
      real(8),intent(in):: v
      real(8):: ans
!
!
!
!        f-function in bremsung function
!
      real*8 g,  t0,  t1, t2, t3
      real*8 t4, t5, t6
      real*8 tmp0, tmp1, dgdv, dgdv2, dgdvs, tmp
!      real*8 cscrn, scrnfp, scrnfb,  smlf1, smlf2 ! don't declare  internal func.

!
!        screening parameter
      g=scrnfb(0,v)
!
      if( g <= 2.d0 ) then
         t0=1.d0-v
         t1=(1.d0+t0*t0)
!        ccz=log(z)/3+f(z); f=coulomb corection function
         t2=(0.25*smlf1(0,g)-ccz)
         t3=(0.25*smlf2(0,g)-ccz)
         tmp0=t1*t2-0.666666666666d0*t0*t3
         if( n == 0 ) then
            ans = tmp0
         else
            dgdv = scrnfb(1,v)
            t4 = .25* smlf1(1,g)
            t5 = .25*smlf2(1,g)
            tmp1 = -2.0*t0*t2+t1*t4*dgdv+0.66666666*(t3-t0*t5*dgdv)
            if( n == 1 )  then
               ans=tmp1
            else
               dgdv2 = scrnfb(2,v)
               dgdvs = dgdv*dgdv*0.25
               ans = 2.0*t2-4.0*t0*t4*dgdv+
     *           t1*(dgdv2*t4+dgdvs*smlf1(2,g) )+1.333333*t5*dgdv-
     *           0.66666666*t0*(dgdv2*t5+dgdvs*smlf2(2,g)  )
            endif   
         endif
!        almost no screening region
      else
         t0=1.-v
         t1=t0*t0-0.6666666*t0+1.
         t2=log(2.*Ebyme*t0/v)-0.5-cscrn(0,g)-fz
         if( n == 0 ) then
            ans = t1*t2
         else
            dgdv=scrnfb(1,v)
            t3=-2.*t0+0.666666
            t4=-1./t0-1./v-dgdv*cscrn(1,g)
            if( n == 1 ) then
               ans = t3*t2+t1*t4
            else
               dgdv2=scrnfb(2,v)
               t5=2.
               t6=1./t0/t0+1./v/v-dgdv *dgdv*cscrn(2,g)-dgdv2*cscrn(1,g)
               ans=t5*t2+2.*t3*t4+t1*t6
            endif
         endif
      endif
      end function bigf
!
!     **********
      function bigg(v) result(ans)
!     **********
!
!     g-function in pair creation function
!
!        screening parameter
!
      real(8),intent(in)::v  
      real(8):: ans

      real(8)::g, tmp
!     real(8):: smlf1, smlf2, cscrn  ! internal func.


      g=scrnfp(v)
      tmp=1.-v
      if( g <= 2.0 ) then 
         ans=(v*v+tmp*tmp)*(.25*smlf1(0,g)-ccz)+.666666*v*tmp*(.25*
     *  smlf2(0, g)-ccz)
      else
         ans=(v*v+tmp*tmp+.666666*v*tmp)*(log(2.*Ebyme*v*tmp)-.5-
     *  cscrn(0,g)-fz   )
      endif
      end function bigg
!     ****************************************************************
!     *                                                              *
!     *  scrnfb:  to give screening parameter for brems              *
!     *  scrnfp:  //                              pair               *
!     *                                                              *
!     ****************************************************************
!
!
!
      function  scrnfb( n, v) result(ans)
      implicit none
#include "Zmass.h"
!
!        gives screening parameter for brems E in gev.
!

      real*8 v
      integer n
      real(8):: ans

      real*8 tmp
      real*8 scrnfp
      tmp=100.0/(Ebyme*Z13)

      select case(n)
         case(0)
            ans=tmp*v/(1.d0-v)
         case(1)
            ans = tmp/(1.d0-v)**2
         case (2)
            ans = tmp*2./(1.-v)**3
         case default
            write(0,*) ' n =',n, ' invalid for scrnfb'
            stop
      end select
      end function scrnfb

      function scrnfp(v) result(ans)
        implicit none
        real(8),intent(in):: v
        real(8):: ans
!
!        gives screening parameter for pair creation
!
        ans=100./(Ebyme*Z13*v*(1.d0-v))
      end function scrnfp
!     ****************************************************************
!     *                                                              *
!     *  smlf1:  auxliary screening function for brems cross-section *
!     *  smlf2:   //                                                 *
!     *  bcorec: correction function for brems at low energies       *
!     *  pcorec: //                      pair  //                    *
!     *  cscrn:  screening function when g>2                         *
!     *                                                              *
!     ****************************************************************
!

!
!
      function smlf1(n,g) result(ans)
      implicit none
!
!        screening function small f1(g)
!
      integer,intent(in):: n  ! 0,1,2
      real(8),intent(in):: g
      real(8):: ans


      real*8  tmp
      real*8 smlf2, bcorec, pcorec, cscrn
      
      if( g < = 1.0d0 ) then
         select case(n)
         case(0)
           ans=((-1.046117*g+2.445063)*g-4.63689)*g+20.83794
         case(1)
           ans = 1.25*g-3.24
         case(2)
           ans =1.25
         case default
           write(0,*) ' n=',n, ' invalid for smlf1'
        end select
      else
         select case(n)
         case(0)
           ans = 19.052795-3.760637*log(g+0.47155)
         case(1)
           ans = -4.184/(g+.952)
         case(2)
           ans = 4.184/(g+.952)**2
         case default
            write(0,*) ' n= ', n, ' invalid in smlf1'
            stop
         end select
      endif
      end function smlf1
!
!
!     ***********
      function smlf2(n,g) result(ans)
      implicit none
      integer,intent(in):: n ! 0, 1, 2
      real(8),intent(in):: g
      real(8):: ans

      if( g <= 1.0d0 ) then
         select case(n)
         case(0)
           ans=((1.49546*g-2.33405)*g-1.73269)*g+20.171278
         case(1)
           ans = -0.172*g-1.930
         case(2) 
           ans = -0.172
         case default
            write(0,*) ' n=',n, 'invalid for smlf2'
         end select
      else
         select case(n)
         case(0)
           ans = 19.052795-3.760637*log(g+0.47155)
         case(1)
           ans = -4.184/(g+.952)
         case(2)
           ans = 4.184/(g+.952)**2
         case default
            write(0,*) ' n= ', n, ' invalid in smlf2'
            stop
         end select
      endif
      end function smlf2
!     ************
      function  bcorec(EE) result(ans)
!     ************
      implicit none
      real(8),intent(in):: EE
      real(8):: ans
!
!        bremsung coreection function at low energy
!         ee=e
!
      ans = 1.+ bcoef/EE
      end function bcorec
!
!     ************
      function pcorec(EE) result(ans)
!     ************
      implicit none
      real(8),intent(in):: EE
      real(8):: ans
!
!         paircreation correction function at low energy
!
      ans = 1.+0.4d-3/(EE-0.99d-3)
      end function pcorec
!     ************
      function cscrn(n,g) result(ans)
!     ************
      implicit none
      integer,intent(in)::n
      real(8),intent(in)::g
      real(8):: ans
!
!        screening function when g.gt.2.
!        approximation made by c(g)=.24/g+.36/g/g
!
      real(8):: tmp
      tmp=1./g
      select case(n)
       case (0)
         ans=g**(-2.2)/(0.20322+2.11737*g**(-1.34562)  )
       case (1) 
         ans = (-2.*.36*tmp-.24)*tmp*tmp
       case(2) 
         ans = (6.*.36*tmp+2.*.24)*tmp*tmp*tmp
       case default
         write(0,*) ' n=',n, 'invalide for cscrn'
         stop
      end select
      end function cscrn
!
!     ************
      subroutine  epBPZpart(Zin)
!     ************
      implicit none
#include "Zglobalc.h"
#include "ZbasicCnst.h"


      real(8),intent(in):: Zin        ! target charge

      real*8 temp
      real*8 epPairLowE

!      real*8 epPair, ccorec  ! internal

      real(8), parameter:: xnorm=0.5d0

      if( Z /= Zin ) then
         Z = Zin
         Z2 = Z*Z
         Z13=Z**(1./3.)
!        coulmb correction function f(z)
         fz=ccorec(Z)
         ccz=log(Z)/3.+fz
         al183z=log(183./Z13)
!        used for bremsung correction at low energies
         bcoef=1.53e-3*sqrt(Z/137.)
!       compute radation length in g/cm2 for virutual matter: A = 1 and z
!       This x0g is used to convert the prob. / X0 into mb
         call epX01(Z, 1.d0, x0g)
      endif
      end subroutine epBPZpart
!     ****************************************************************
!     *                                                              *
!     *  ccorec: coulomb correction function f(z)                    *
!     *                                                              *
!     ****************************************************************
!
!
!        coulomb correction function  a=z/137 ccorec=a**2 sigma(1/((n**
!        +a**2) from n=1 to inf .  approx formula used.
!
!
      function ccorec(z) result(ans)
      implicit none
      real(8),intent(in):: z
      real(8):: ans

      real*8 a
!

      a=z/137.d0
      a=a*a
      ans = ( ( (-0.002*a+0.0083)*a-0.0369)*a+0.20206+1./(1.+a)   )*a
      end function ccorec

!      end module BPKochMotz
      end module BPPS

