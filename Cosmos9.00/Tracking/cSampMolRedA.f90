! This is a copy from Epics.  Name is changed with C instead of ep
! In future the one in Epics should be removed and should use this.
!
!old  test program is located prog/Elemag/Util/Moliere/testSampMolRedA.f90
module cmodSampMolReducedA
  implicit  none


  real(8),parameter::xmax=200.d0  ! max x to be sampled
      ! sampling function f= f0 + f1/B + f2/B^2.  f0=2exp(-x)
      ! f1 and f2 become negative in some region.  
      ! The boundries set to the point where f1 or f2 change the
      ! sign.    
  real(8),parameter::Txmax=50.d0  ! below this, f1,f2 use  tables
  integer,parameter:: nregion=7   ! # of such regions
  integer,parameter:: nregionA=nregion + 1   ! # of  regions upto xmax
      !  there is a very small gap between two consecutive regions.
  complex(8),save:: xr(nregionA) = &   ! +1 : last one is by formula
         ! f1>0                    f1<0
         ! f2>0                    f2>0  
       (/(0.d0, 0.31517d0), (0.31528d0,  0.35880d0), &
         ! f1<0                 f1<0                   f1>0
         ! f2<0                 f2>0                   f2>0
       (0.35892d0, 1.8660d0), (1.8662d0, 3.2440d0), (3.2443d0, 5.0400d0),  &
         !  f1>0                   f1>0
         !  f2<0                   f2>0
       (5.0405d0,  14.398d0), (14.399d0,  Txmax),  (Txmax,  xmax)/)   
!  at x>Txmax f1,f2= a(x/Txmax)**b;  a for f1 is f1a etc
  real(8),parameter:: f1a = 0.000868d0 
  real(8),parameter:: f1b = -2.0598d0  ! fit error  mostly < 0.2% max 1
  real(8),parameter:: f2a =  0.0001508d0 
  real(8),parameter:: f2b = -2.6645d0  ! mostly <0.5 (x<100) % 0.5~2.3%  (x>100)   This is not important.

  integer :: region   ! selected region index
  real(8):: x1, x2  ! lower and upper boundary of a region

 !    integral value of f0=2exp(-x), f1 and, f2 in each region
  real(8),save:: FI(0:2,nregionA)
  real(8):: Ft(nregionA) ! total integral vlaue F0(i) + F1(i)/B + F2(i)/B/B 
  real(8):: Ftall   ! sum of Ft
  real(8):: Ftpos   ! sum of   positive f0, f1/B f2/B/B  in a given region
  integer:: negaflag  ! # of negative  FI(i,j) for a given (i,j)
  integer:: funcno    ! selectef function number 0, 1, or 2  for sampling
  real(8):: add   ! working var. to  sum up quantities
  real(8):: B2   ! B^2
  real(8):: Pos, Neg, Tot
  real(8):: eff  ! used for sampling  when some of f1,f2 becoem negative
                 ! see later
!     to store (x, fi) for number of x's. (i=1,2)
  type xy
     real(8), allocatable::xa(:), ya(:)
     integer::n
     integer::id   ! id for sampling
  end type xy
  type(xy),save:: f1(nregion), f2(nregion)  ! table of f1,f2 in each region

!           for f1, use different xbin in the following region
  integer,parameter:: nbin1=7
  real(8),save::xbinr1(nbin1)=(/0.25d0, 0.33d0, 1.0d0, 2.2d0, 4.5d0, 5.1d0, xmax/)
        ! for increasing x,  xbin1(i) is used if x is below xbinr1(i)
!  real(8),save::xbin1(nbin1)=(/0.01d0, 0.001d0, 0.01d0, 0.05d0, 0.01d0, 0.025d0, 0.5d0/) 
  real(8),save::xbin1(nbin1)=(/0.02d0, 0.002d0, 0.02d0, 0.05d0, 0.1d0, 0.05d0, 0.5d0/) 
!           same for f2
  integer,parameter:: nbin2=8
  real(8),save::xbinr2(nbin2)=(/0.25d0, 0.33d0, 2.5d0, 4.1d0, 7.0d0, 13.0d0, 17.0d0,  xmax/)
            ! for increasing x,  xbin2(i) is used if x is below xbinr2(i)
!  real(8),save::xbin2(nbin2)=(/0.01d0, 0.001d0, 0.01d0, 0.05d0, 0.01d0, 0.05d0,    0.01d0,  0.5d0/) 
  real(8),save::xbin2(nbin2)=(/0.02d0, 0.001d0, 0.02d0, 0.05d0, 0.05d0, 0.1d0,    0.1d0,  0.5d0/) 


  contains

  subroutine cSampMolReducedA(B, thetasq, icon)
    use modcsampAF
    implicit none
    real(8),intent(in)::B  ! Moliere's B 
    real(8),intent(out)::thetasq ! sample reduced angle square (x)
               ! tetasq = x ; teta= sqrt(tetasq*B*xc2)
    integer,intent(out)::icon  ! =0, ok, =1, not ok-
    real(8):: u, ans
    real(8)::x   !  thetasq
    logical,save::first=.true.
    integer i
!    real(8):: cf1approx, cf2approx  



    if(first) then
       call CSampMolini
       first = .false.
    endif
!   sampling:  first fix the region from which x is sampled
!
    B2 = B*B
    do i = 1, nregionA
       Ft(i) = FI(0, i) + FI(1,i)/B + FI(2,i)/B2
    enddo
    Ftall = sum(Ft(:))
!///////////
!    write(0,*) 'Moliere func B=',B, ' Ftotal =', Ftall 
!/////////////
    call rndc(u)
    add = 0.
    do i = 1, nregionA
       add = add + Ft(i)
       if( u <=  add/Ftall )  exit
    enddo
    region = i   ! selected region
!
!     fix which function f0, f1, f2


    Ftpos = 0.
    do i = 0, 2
       if( FI(i, region) > 0.) then
          Ftpos = Ftpos  + FI(i, region)/B**i
       endif
    enddo

    negaflag = 0
    do i = 1, 2
       if( FI(i, region) < 0.) then
          negaflag = negaflag + 1
       endif
    enddo

    do while (.true.)
       call rndc(u)
       add = 0.
       do i = 0, 2
          if( FI(i, region) >= 0.) then
             add = add + FI(i, region)/B**i
             if( u <= add/Ftpos) then
                funcno = i
                exit
             endif
          endif
       enddo
!    funcno function is used 
       if(funcno == 0 ) then
        ! 2exp(-x) dx  in region.   solv
        ! (exp(-x1) -exp(-x2))u = (exp(-x1) - exp(-x))
          x1 = real( xr(region))
          x2 = aimag( xr(region ) )
          call rndc(u)
          x = -  log(exp(-x1) - (exp(-x1)- exp(-x2))*u)
       elseif(funcno == 1) then
          call csampf1(x)
       else
          call csampf2(x)
       endif
       if(negaflag > 0 ) then
        !  negative FI(1:2,region ) exist
        !  rejection using eff below

        !   S + N = T  (S; sum of func's  used for sampling
        !               N: sum of  negative func.
        !               T: true function       
        !  We demand  eff*S = T ; i.e eff*S =S + N (eff=efficiency)
        !               eff = 1 + N/S < 1.0 
        ! 
        !  eff*(sum FI>0) = FI(0,region) + FI(1,region)/B + FI(2,region)/B2
        !             = Ft(region)
        !  eff = Ft(region)/(sum FI>0)
        !  accept if (u< eff) 
          Pos = 0.  ! S
          Neg = 0.  ! N
          Tot = 0.  ! T
          Pos = 2*exp(-x) 
          call cf1approx(x, ans)
!           call  csampAFIntp(f1(region)%id, x, ans)
          if( FI(1, region) > 0. ) then
             Pos = Pos + ans/B
          else
             Neg = Neg + ans/B
          endif
          
!          call  csampAFIntp(f2(region)%id, x, ans)
!          ans = cf2approx(x)
          call cf2approx(x, ans)
          if( FI(2, region) > 0.) then
             Pos = Pos + ans/B2
          else
             Neg = Neg + ans/B2
          endif
          Tot = Pos + Neg
          eff = Tot / Pos
          call rndc(u)
          if( u < eff ) exit
       else
          exit
       endif
!        We must repeat until accepted.
    enddo
    icon = 0
    thetasq = x
  end subroutine cSampMolReducedA

  subroutine cSampMolini
    use modcsampAF
    implicit none


    real(8),external:: cMolfunc1, cScotD2


    real(8):: error
    real(8),save:: eps=1.d-5
    integer::jcon
    real(8)::x   !  thetasq
    real(8)::xlast
    logical,save::first=.true.
    integer:: i, j, k, n
    real(8):: temp
!      since this is one pass routine, if called 
!     > 1, do nothing.
    if(.not. first) then
       return
    else
       first = .false.
    endif
!         integrate f0,f1,f2 in each region


    do i = 1, nregionA  ! entire region
       x1 = real( xr(i) )
       x2 = aimag( xr(i) )
!           f0 = 2exp(-x); integral = -2exp(-x)
       FI(0, i) =  2.d0*(exp(-x1) - exp(-x2)) 
!       write(0,*) 'region ', i, ' FI(0,i)=',FI(0, i)
       if( i < nregionA ) then
          call kdexpIntF(cMolfunc1, x1, x2, eps, FI(1,i), error, jcon)
       else
          !  a/(b+1) (x/x1)**(b+1) from x1,x2  
          temp = f1b + 1.0d0   !< 0.
          FI(1,i) =x1* f1a/temp*( (x2/x1)**temp - 1.d0)
       endif
!       write(0,*) 'region ', i, ' jcon =', jcon, ' error=', error, ' FI(1,i)=',FI(1,i)
       if(i < nregionA) then
          call kdexpIntF(cScotD2, x1, x2, eps, FI(2,i), error, jcon)
       else
          temp = f2b + 1.0d0  
          FI(2,i) = x1*f2a/temp*( (x2/x1)**temp - 1.d0)
       endif

!       write(0,*) 'region ', i, ' jcon =', jcon, ' error=', error, ' FI(2,i)=',FI(2,i)
    enddo
!          compute f1 in each region; count number  of needed points
!          in each region and then allocate size
    k = 1   ! bin size index  

    do i = 1, nregion
       x = real( xr(i) )
       x2 =aimag( xr(i) )
       n = 0   ! # counter 
       do while (x <= x2)
          if(x > xbinr1(k) ) then
             k = k + 1
             if( x >  xbinr1(k) ) then
                write(0,*) ' strange'
                stop
             endif
          endif
          n = n + 1
          xlast = x
          x = x + xbin1(k)
       enddo
       if( xlast < x2-xbin1(k)*0.1d0 ) then
          n = n + 1
       endif
!///////////
!       write(0,*) ' allocate func1 region ',i, ' n=', n
!/////////////////
       allocate( f1(i)%xa(n), f1(i)%ya(n) )
       f1(i)%n = n
    enddo
!           put data in each region
    k = 1

    do i = 1, nregion
       x = real( xr(i) )
       x2 =aimag( xr(i) )
       n = 0
       do while (x <= x2)
          if(x > xbinr1(k) ) then
             k = k + 1
          endif
          n = n + 1
          f1(i)%xa(n) = x
          f1(i)%ya(n) = cMolfunc1(x)
          xlast = x
!////////////////////
!          write(*,'(a, 2i4, 1p,2g13.4)') &
!               'f1  reegion ', i, n,  f1(i)%xa(n), f1(i)%ya(n)    ! only at test time
!////////////////////
          x = x + xbin1(k)
       enddo
       if( xlast < x2-xbin1(k)*0.1d0 ) then
          n = n + 1
          f1(i)%xa(n) = x2
          f1(i)%ya(n) = cMolfunc1(x2)
!          write(0,'(a, 2i4, 1p,2g13.4)') &
!               'f1  region ', i,n, f1(i)%xa(n), f1(i)%ya(n)    ! only at test time
       endif
    enddo
!        same for f2

!          compute f2 in each region; count number  of needed points
!          in each region and then allocate size

    k = 1

    do i = 1, nregion
       x = real( xr(i) )
       x2 = aimag( xr(i) )
       n = 0
       do while (x <= x2)
          if( x > xbinr2(k) ) then
             k = k + 1
             if( x >  xbinr2(k) ) then
                write(0,*) ' strange'
                stop
             endif
          endif
          n = n + 1
          xlast = x
          x = x + xbin2(k)
       enddo
       if( xlast < x2-xbin2(k)*0.1d0 ) then
          n = n + 1
       endif
!///////////
!       write(0,*) ' allocate func2 region ',i, ' n=', n
!/////////////////
       allocate( f2(i)%xa(n),  f2(i)%ya(n) )
       f2(i)%n = n
    enddo
!           put data in each region
    k = 1

    do i = 1, nregion
       x = real( xr(i))
       x2 = aimag( xr(i))
       
       n = 0
       do while (x <= x2)
          if(x > xbinr2(k) ) then
             k = k + 1
          endif
          n = n + 1
          f2(i)%xa(n) = x
          f2(i)%ya(n) = cScotD2(x)
!////////////////////
!          write(*,'(a, 2i4, 1p,2g13.5)') &
!               'f2  region ', i,n, f2(i)%xa(n), f2(i)%ya(n)    ! only at test time
!////////////////////
          xlast = x
          x = x + xbin2(k)
       enddo
       if( xlast < x2-xbin2(k)*0.1d0 ) then
          n = n + 1
          f2(i)%xa(n) = x2
          f2(i)%ya(n) = cScotD2(x2)
!////////////////////
!          write(*,'(a, 2i4, 1p,2g13.5)') &
!               'f2  region ', i,n, f2(i)%xa(n), f2(i)%ya(n)    ! only at test time
!////////////////////
       endif

    enddo
       !           initialzize sampling for f1,f2 in each  region if F >0
    do i = 1, nregion
       if(FI(1,i) > 0. ) then
          call csampAF0byArray(f1(i)%xa, f1(i)%ya, f1(i)%n, f1(i)%id)
       else
!              only for getting function value
          call csampAF0byArray_b(f1(i)%xa, f1(i)%ya, f1(i)%n, f1(i)%id)
       endif
       if(FI(2,i) > 0. ) then
          call csampAF0byArray(f2(i)%xa, f2(i)%ya, f2(i)%n, f2(i)%id)
       else
          call csampAF0byArray_b(f2(i)%xa, f2(i)%ya, f2(i)%n, f2(i)%id)
       endif
    enddo
  end subroutine cSampMolini

!  function epf1approx(x)
  subroutine cf1approx(x,ans)
    use modcsampAF
    implicit none
!    real(8):: cf1approx
    real(8),intent(in):: x    ! reduced angle sq.

    real(8)::ans
    if(region < nregionA) then
       call  csampAFIntp(f1(region)%id, x, ans)       
    else
       ans = f1a*(x/Txmax)**f1b
    endif
!    epf1approx = ans
  end subroutine cf1approx

!  function epf2approx(x)
  subroutine cf2approx(x,ans)
    use modcsampAF
    implicit none
!    real(8):: epf2approx
    real(8),intent(in):: x    ! reduced angle sq.

    real(8)::ans
    if(region < nregionA) then
       call  csampAFIntp(f2(region)%id, x, ans)       
    else
       ans = f2a*(x/Txmax)**f2b
    endif
!    epf2approx = ans
  end subroutine cf2approx


  subroutine csampf1(x)
    use modcsampAF
    implicit none
    real(8),intent(out):: x    !sampled reduced angle sq.

    real(8):: u, temp
    if(region < nregionA) then
       call csampAF(f1(region)%id, x)
    else
       !  a/(b+1)(x/x1)**(b+1) ; 
       !  solve:  ( (x2/x1)**(b+1)  - 1.) u = (x/x1)**(b+1) -1
       !  ( (x2/x1)**(b+1)  - 1.) u +1. = (x/x1)**(b+1) 
       !  x = x1(1. +[ (x2/x1)**(b+1) -1]u)**(1./(b+1))
       call rndc(u)
       temp = f1b+1.0d0 
       x= Txmax*( ( (xmax/Txmax)**temp -1.d0)*u + 1.d0 )**(1.d0/temp)
    endif
  end subroutine csampf1

  subroutine csampf2(x)
    use modcsampAF
    implicit none
    real(8),intent(out):: x    !sampled reduced angle sq.
    real(8):: u, temp
    if(region < nregionA) then
       call csampAF(f2(region)%id, x)
    else
       call rndc(u)
       temp = f2b+1.0d0 
       x= Txmax* ( ( (xmax/Txmax)**temp -1.d0)*u + 1.d0 )**(1.d0/temp)
    endif
  end subroutine csampf2

end module cmodSampMolReducedA
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

  
