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
