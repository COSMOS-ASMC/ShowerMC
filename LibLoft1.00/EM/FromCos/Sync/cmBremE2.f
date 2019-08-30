!       test cmBremE2: second term of mag. brem.
!	  implicit none
!	  real*8  v, up
!	  integer j, i
!
!          read(*,*) up
!
!          do  i = 1, 50000
!             call cmBremE2(up, v, j)
!             write(*, *) sngl(v), j
!          enddo
!	  end
!
!	     samples fractional gamma ray energy from 
!       v phai(2zeta)dv which is the 2nd term of the
!       syncroton emission, where phai(x) = x K2/3(x).
!       This can be rewritten as
!      18Uz
!   ------------ z K2/3(z). Here U is upsilon.
!    (2+3Uz)^3
!
!    z^2 K2/3(z)  and (2+3Uz)^(-3) are approximated by
!    a number of power functions.
!    
!   z^2 K2/3(x)  by 5 
!               nOte:  --> 1.0747764 z^(4/3)   z-->0
!                      -->  sqrt(pi/2) z^(3/2) exp(-z)   z-->inf.
!    The latter func. has peak value, 0.513 at z=3/2
!    and the slope at z is 3/2 - z.  
!
!   1/(2+3Uz)^3 --->  1/8.  z-->0
!                     1/z^1.5     with a value of 1/64 at z = 2/(3U).
!                     1/(3Uz)^3  z--> inf.
!
!   
      subroutine cmBremE2(up, v, nc)
  	  implicit none
	  real*8 up       ! input Upsilon value
	  real*8 v        ! output. sampled fractional value
          integer nc      ! output. number of rejections. average trial
!                                   number is 1.16 for up=0.01 and 1.26 for
!                                   up=10.
!
	  integer n1, n2, n, i, na
	  parameter (n1 = 3, n2 = 5, n= 12)
	  real*8 c1(n1), x10(n1), p1(n1), x1p(n1+1), coef1(n1)  ! for 1/(2+3Uz)^3
!
	  real*8 c2(n2), x20(n2), p2(n2), x2p(n2+1), coef2(n2) ! for z^2K2/3(z)
!
	  real*8 coef(n), p(n), xp(n+1)
	  real*8 z, fz, ft, ck23, u
	  logical first

	  data p1/0.04433d0, 1.5d0, 2.8125d0/
	  data c1/0.1195d0,  1.5625d-2, 3.05175d-5/

!
	  data x20/0.01d0, 1.d0, 2.511886d0, 7.4131d0,
     *        19.95d0/
	  data c2/2.3096d-3, 4.9463d-1, 0.41808d0,
     *        1.5448d-2, 2.427d-7/
 	  data p2/-1.3298d0, -0.446d0, 1.04103, 5.92368d0,
     *        18.007d0/
     

	  data first/.true./
	  save first, x1p, coef1, x2p, coef2, x10, x20
	  
!

	  x10(1) = 0.01d0/up
	  x10(2) = 2.d0/3.d0/up
	  x10(3) = 10.d0/up
!

	  call ksampPwX(c1, p1, x10, n1, x1p(2), coef1 )
	  x1p(1) = 0.
	  x1p(n1+1) =  30./up
!          write(*,*)' x1p, coef1=', x1p, coef1
!	
	  if(first) then
		 call ksampPwX(c2, p2, x20, n2, x2p(2), coef2 )
  		 first = .false.
	     x2p(1) = 0.
	     x2p(n2+1) = 30.
!             write(*,*) ' x2p, coef2=', x2p, coef2
      endif	
	call ksampPwMrg(coef1, p1, x1p, n1,
     *  coef2, p2, x2p, n2,
     *  coef,  p,  xp,  na)
!       write(*,*)'n, coef, p, xp=',  na, coef, p, xp
	do i =1, 1000
	   call ksampPw(i, coef, p, xp, na, z, fz)
!           write(*,*) ' z, fz=', z, fz
	   ft = z*z * ck23(z) / (2.+ 3.* up*z)**3
	   call rndc(u)
	   if(u .lt. ft/fz) goto 10
	enddo
 10	continue	
	nc = i
	v =  3.*up*z/(2.+ 3.*up*z)
	end

