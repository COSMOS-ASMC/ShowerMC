!  This is copy of epMolfunc.f90 and epMolBlogB.f90
!  in Epics.  epxxxx is changed to cxxxx
!  In future, Epics should call those here and remove epMolfunc.f90 
!            and epMolBlogB
!  kpolyg and   keimlg are also Epics functions
!  so they are copied to Cosmos/KKlib and renamed kpolygC and keimlgC
!
!   ifort  cMolfunc.f90  -L$EPICSTOP/lib/PCLinuxIFC64 -lepics
!   -L$COSMOSTOP/lib/PCLinuxIFC64 -lcosmos

!    compute Moliere's F0, F1, F2 
!  which is
! 2pi  f(theta) theta dtheta = ( F0(x)  + F1(x)/B   F2(x)/B^2) dx
!                 F0= 2exp(-x)
!                 F1 = cMoliere1(x)
!                 F2 = cScotD2(x) (don't use cMoliere2)
!   x = teta^2 (teta reduced angle; rad^2)
!  
!
!program main
!implicit none
!real(8):: teta, x
!real(8):: cMolfunc1, cScotD2
!real(8):: f0, f1, f2, B
!!
!!/////////////
!B=1.
!!////////////
!do while (B< 30.)
!   teta = 0.
!   do while( teta<14.)
!      x = teta**2
!      f0 = 2*exp(-x)
!      f1 = cMolfunc1(x)
!      f2 = cScotD2(x)
!      write(*,'(1p,7g13.4)') B, teta, x, f0, f1/B, f2/B**2, f0+f1/B+f2/B/B
!      teta = teta + 0.001d0
!   enddo
!!//////////////
!   if(B == 1.0) exit
!!/////////////////
!   write(*,*) " "
!   B = B + 2.
!enddo
!end program main
!! program main
!   write(*,*) teta, cMolfunc1(x), x
!   teta = teta + 0.1d0
!enddo
!end

! ifort cMolfunc.f90 -L$EPICSTOP/lib/PCLinuxIFC64 -lepics
!   cMolfunc2 based on Bethe's paper dose not give 
! the correct result so we should use cScotD2 instead.
!
!!  test cMolfunc1

!implicit none
!real(8):: teta, x
!real(8):: cMolfunc1
!
!teta = 0.
!do while( teta<11.)
!   x = teta**2
!   write(*,*) teta, cMolfunc1(x), x
!   teta = teta + 0.1d0
!enddo
!end

!!  test cMolfunc2
!  program main
!   implicit none
!   real(8):: teta, x
!   real(8):: cMolfunc2
!   
!   teta = 0.
!   do while( teta<3.7)
!      x = teta**2
!      write(*,*) teta, cMolfunc2(x), x
!      teta = teta + 0.1d0
!   enddo
! end program main
!!  test cMolfunc2I and  cMolfunc2I2
! program main
!  implicit none
!  real(8):: teta, x
!  real(8):: cMolfunc2I
!  real(8):: cMolfunc2I2
!  real(8):: ans1, ans2 
!  teta = 0.
!  do while( teta<1.)
!     x = teta**2
!     ans1 = cMolfunc2I(x)
!     ans2 = cMolfunc2I2(x)
!      write(*,*) teta, ans1, ans2, x
!     teta = teta + 0.1d0
!  enddo
! end program main
!!
function cPsi(n, x)
  ! kpolyg is  n-th derivative of d log(G(x))/dx ( G(x) is the
  ! gamma function.
  ! On the other hand, 
  ! Bethe's paper: Phys. Rev. 59, 1953
  ! uses  psi(x) =  d log(G(x+1))/dx.  and its derivative
  ! To avoid confusion, we use cPsi which calls kpolyg
  !   cPsi(n, x) = kpolyg(n, x+1)
  ! note: 
  ! cPsi(0, x) = 1/x + kpolyg(0, x) = kpolyg(0,x+1)  
  implicit none
  real(8)::cPsi
  integer,intent(in):: n  ! 0,1,2...  n-th derivative
  real(8),intent(in):: x  !

  real(8),external:: kpolygC
  cPsi = kpolygC(n, x+1.0) 
end function cPsi


  

function cMolfunc1(x)
  ! Moliere's f1(x) as defined by Eq.(28) of
  ! Bethe's paper: Phys. Rev. 59, 1953
  !*** Tested ***
  implicit none
  real(8):: cMolfunc1
    !   = 2exp(-x)(x-1)(Ei(x) - ln(x)) - 2(1-2exp(-x))
  real(8),intent(in)::x   ! reduced angle 
  real(8),external::keimlgC   !  Ei(x) - ln(x)

  real(8):: temp

  temp = 2* exp(-x)
  cMolfunc1 = temp*(x-1.)* keimlgC(x) - 2*(1.0-temp)
end function cMolfunc1

! function cMolfunc2(x)
! implicit none
!   ! Moliere's f2(x) as defined by Eq.(29) of
!   ! Bethe's paper: Phys. Rev. 59, 1953
!   !
!   real(8):: cMolfunc2
!     !   = {psi(2)^2 + psi'(2)}(x^2-4x+2)
!     !           + integral()
!     !  integrand is complex and the integration can be
!     !  expanded as infinite series. 
!   real(8),intent(in)::x   ! reduced angle 
!   real(8)::cMolfunc2I
!   real(8)::cMolfunc2I2
!   logical,save:: first = .true.
!   real(8),save:: psi2sq
!   real(8),save:: psip2
!   real(8):: cPsi
! 
!   if(first) then
!      psi2sq = cPsi(0, 2.d0)**2
!      psip2 =  cPsi(1, 2.d0)
!      first = .false.
!   endif
! 
!   cMolfunc2 = (psi2sq + psip2)* (x*(x-4.)+2.0) +  cMolfunc2I2(x) 
! 
!  
! !  cMolfunc2 =( (psi2sq + psip2)* (x*(x-4.)+2.0) +  &     cMolfunc2I(x) ) 
! 
! end function cMolfunc2
! 
! function cMolfunc2I(x)
! implicit none
!   !  func2 includes an integral which can be expressed 
! !    by a series which is this funcion.
!   real(8):: cMolfunc2I
!   real(8),intent(in):: x
!   
!   real(8):: cPsi
! 
!   integer::n
!   real(8),parameter:: eps= 1.e-5
!   real(8):: sum, fac1, fac2, fac3,  temp, nd
!   real(8):: term, x1, x2, x3, psi2, term2
!   real(8),save:: Euler = 0.57721566d0
!   
!   logical,save::first=.true.
!   real(8),save::fact(0:50)
!   integer i
!   if( first ) then
!      fact(0) =1.
!      do i = 1, 50
!         fact(i) = fact(i-1)*i
!      enddo
!      first = .false.
!   endif
!   
!   sum = 0.
!   psi2 = cPsi(0, 2.d0)
!   do n = 0, 30
!      temp = n+1
!      nd = n
!      term2 = x**(n+1)/fact(n+1) * ( x/(n+2)*(x/(n+3) -2.0) + 1.)
!      term = (cPsi(0, nd) + Euler - psi2)/temp  * term2
!      sum = sum + term
!      if( term == 0.) exit
!      if( abs(term/sum) < eps) exit
!   enddo
!   cMolfunc2I = sum
! end function cMolfunc2I
! 
! function cMolfunc2I2(x)
! implicit none
!   !  func2 includes an integral which can be expressed 
! !    by a series which is this funcion.
!   real(8),intent(in):: x
! 
!   real(8):: cMolfunc2I2
! 
! 
!   real(8):: cMolfunc2I2f
!   external  cMolfunc2I2f
!   real(8),save:: eps = 1.d-4
!   real(8):: ans, error
!   integer icon
! 
!   real(8):: xx
!   common /temp/ xx
! 
!   xx =x
!   call kdexpIntF(cMolfunc2I2f, 0.d0, 1.d0, eps, ans, error, icon)
!   
!   cMolfunc2I2 = ans
! end function cMolfunc2I2
! 
! function cMolfunc2I2f(tt)
! implicit none
! 
!   real(8):: cMolfunc2I2f
!   real(8),intent(in):: tt(2)
!   
!   real(8):: eps, psi2, sing, sing2, t2, t
!   real(8):: cPsi
!   real(8),parameter::small = 1.d-4
!   real(8):: xx
!   common /temp/ xx
!   real(8):: x
!   
!   x = xx
!   psi2 = cPsi(0, 2.d0)
!   t2 = tt(2)
!   if( t2 < 0. ) then
! !     t = 0.0- t2  !  1-t = 1+t2
!      t = -t2  
!      sing = (log(t)/(1.+t2) -psi2 ) 
!      if( t < small) then
!         sing2 = (x-x**2 + x**3/6.)
!      else
!         sing2 =  &
!          ( (1.0+t2)**2 * exp(x*t) - 1. -(x-2)*t -(x**2/2 -2*x +1)*t**2 ) /t**3
!      endif
!      cMolfunc2I2f = sing* sing2
!   elseif( t2 >= 0.) then
!      t = 1.-t2
! 
!      if( t2< small) then
! !            log(1-t2)/t2 = -1
!         sing =  -1.- psi2
!      else
!         sing =  (log(t)/t2 -psi2 )
!      endif
!      cMolfunc2I2f = sing * &
!           ( t2**2 * exp(x*t) - 1. -(x-2)*t -(x**2/2 -2*x +1)*t**2)
!   endif
! end function cMolfunc2I2f
! 
! ifort cScotD2.f90 cMolfunc.f90  -L$EPICSTOP/lib/PCLinuxIFC64 -lepics -L$COSMOSTOP/lib/PCLinuxIFC64 -lcosmos
!program main
!implicit none
!real(8)::x, teta, f2 
!real(8)::cScotD2
!teta = 0.0d0
!do while (teta< 13.d0)
!   x = teta ** 2
!   f2 = cScotD2(x)
!   write(*,*) teta, f2, x
!   teta = teta + 0.1d0
!enddo
!end program main
!
!
function cScotD2(x)
  implicit none
! Moliere's F2(x) :SCOT Rev.Mod.Phys. (A24)
!  (p.311). with alpha=3, beta=1. 
! 
!       at around x=teta^2 =4.2^2, the accuracy
!   is  worst ( 3 digit accuracy).
!   at others, much more accurate.
!   *** Tested
! With Bethes forumulat in Phys. Rev. I couldn't
! get correct result although numerical table
! by Bethe and Scot coinside.
!
!  G(1) SUM  G(3+k)(-x)^k/(k!G(1+k))
!      x   (  Psi^2(0,2+k) + Psi(1,2+k)  )
!  SUM is over k=0 to inf.
!  G(x) = gamma(x)
!  Psi(n,x)  = kpolyg(n,x+1)    
  real(8)::cScotD2
  real(8):: x  ! reduced theta^2
  
  real(8),external::kpolygC, kgamma
  real(8),parameter:: bigx=4.21*4.21
  integer k

  logical,save:: first= .true.
  real(8):: fac, psi1, psi2, xk, sum, term, termold
  real(8):: temp, temp2, temp3, xbyfac, facbyx, temp1

  real(8):: cPsi
  real(8),parameter:: eps= 1.d-5
  
  
  if( x < bigx ) then
     
!  G(1) SUM  G(3+k)(-x)^k/(k!G(1+k))
!      x   (  Psi^2(0,2+k) + Psi(1,2+k)  )

!  SUM k=0, 1, ..
!    (k+2)! /k!/k! = (k+2)(k+1)/k! 
     psi1 = cPsi(0, 2.d0)
     psi2 = cPsi(1, 2.d0)
     sum = 2.* (psi1**2 + psi2)
     xbyfac = 1.
     if(x /= 0. ) then
        do k = 1, 200
!           k+2 :  k+1
!           psi2 = psi2 - 1.d0/(k+1)**2
           temp = k + 3
           temp2 = k + 1
           temp3 = k+2
           psi1 = psi1 + 1.d0/temp3
              !           psi1 = cPsi(0, temp3)
              !           psi2 = cPsi(1,temp3)
           psi2 = psi2 - 1.d0/temp3**2
              !           term = kgamma(temp)* (-x)**k/ kgamma(temp2)**2  * &
           xbyfac = xbyfac*(-x)/k
           !           term =temp3*temp2 * (-x)**k/ kgamma(temp2)  * &
           term =temp3*temp2 * xbyfac  *   (psi1**2 + psi2)
           sum = sum + term
           if( abs(term /sum) < eps) exit
        enddo
     endif
  else
! For large x; m=2, beta=1; use asymptotic expansion
!   D2(b+m, b, -x) = 2(-1)G(1)/x^3 SUM
!           x [ G(3+k)(2+k)!/(k! x^k  )]
!           x [ psi(0, 2+k) + psi(0,2+k)-log(x) ]
!   psi(k+2) = psi(k+1)+1/(k+1)
!  SUM from k=0 to n; n some artibtray number.
!   k =0;
     psi1 = cPsi(0, 2.d0)  ! k+2

     sum = 4.0 * (2* psi1  - log(x))
     facbyx = 2.
     do k = 1, 50
!        facbyx = (k+1)*facbyx/x
        temp1 = k + 1
        temp2 = k + 2
        temp3 = k + 3
!        facbyx = kgamma(temp3)/x**k
        facbyx = facbyx*temp2/x
        psi1 = cPsi(0, temp2)
        term = facbyx * (k+2)*(k+1)*(2*psi1  - log(x))
        if(k>8) then
           if(abs(term) > abs(termold)) exit
        endif
        termold = term
        sum = sum + term
        if( abs( term/ sum) < eps )  exit
     enddo
     sum = -2* sum /x**3
  end if
  cScotD2 = sum
end function cScotD2
  
  
  
!program main
!implicit none
!real(8):: sb  
!real(8):: B  
! sb = 1.1d0
! do while (sb < 35.)
!    call cMolBlogB(sb, B)
!    write(*,*) sb, B, (B-log(B)-sb)/B
!    sb = sb+0.05d0
! enddo
!end program main
!

subroutine cMolBlogB(sb, B)
implicit none
real(8),intent(in):: sb  ! solve  B- log(B) = sb for given sb
real(8),intent(out):: B  ! 

real(8)::sbb, fB, fBp
if(sb <= 1.d0) then
   write(0,*) ' b for B-log(B) = b is too small b=',b
   stop
endif
if(sb < 5.d0) then
!     sb<2.5-->error <0.06 %
!       <5     error <0.005% 
   sbb =sb
   B = -0.591837 + (2.06538 -0.11690*sbb)*sbb
   fB = B - log(B) -sb
   fBp = 1.d0 - 1.d0/B
   B =B-  fB/fBp
!   sbb = B - log(B)
!   B = -0.591837 + (2.06538 -0.11690*sbb)*sbb
elseif( sb < 15.d0) then
!  5<sb<15;  next one;  |error|  mostly 0.1% max 0.3%
   B = 1.154804 + (1.183689 -0.00462208*sb)*sb
!    by next, error becomes  almost 0
   fB = B - log(B) -sb
   fBp = 1.d0 - 1.d0/B
   B =B-  fB/fBp
elseif( sb < 50.d0 ) then
!   |error| < 0.02 % max 0.05 %
   B = 1.9368162469 +(1.074931428 -0.0007441876*sb)*sb
!    by next almost 0
   fB = B - log(B) -sb
   fBp = 1.d0 - 1.d0/B
   B =B-  fB/fBp
else
   write(0,*) ' sb for B-log(B)=sb too large=',sb
   stop
endif
end subroutine cMolBlogB

  
