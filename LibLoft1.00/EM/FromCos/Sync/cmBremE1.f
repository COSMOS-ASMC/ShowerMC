      subroutine cmBremE1(up, x)
!
!     Samples synchrotron gamma energy as the fraction
!     to the electron energy.  Test for up = 10**-8 to 100.
!     This is for the first part of Eq.5 Brainerd and Petrosian
!     (APJ 320. 1987)
      implicit none
!
      real*8 up !  input. Upsilon value. = Ee/Me *  B/Bcr
!
      real*8 x  !  output.  sampled fractional gamma energy

!

      integer n1, n2, n,  i, nx
      parameter (n1 = 4, n2=3, n=8) 
      real*8  c(n), p(n), xp(n+1)
      real*8  fx,  ft
      real*8  x01(n1),  x02(n2), c1(n1), c2(n2)
      real*8  p1(n1), p2(n2)
      real*8  xp1(n1+1), coef1(n1)
      real*8  xp2(n2+1), coef2(n2)
      real*8  z, u, cmBremF11
      data c1/2.149, 0.918, .5119, 0.00877/
      data p1/-0.3333333333d0, 0.0d0, 1.003, 5.54/
      data x01/1., 1.0, 1.333d0, 5.9566/
      data xp1/0.0, 7.795048057547918E-02, .7446141661381438,
     *        3.3841299222059, 15.0/
      data coef1/2.149, .918, .6829513524709896, 172.3816110479136/


      data p2/1., 1.666666, 3./
      data c2/1., 0., 0./
      data x02 /1., 0., 1./
      
!      --------------------------------------------------
      c2(2) = up/3.
      c2(3) = 1./9/up/up

      x02(2) = 1./3./up

!       ******************* next is no more needed **********
!       since the result is put in the data statement.
!
!      call ksampPwX( c1, p1, x01, n1, xp1(2), coef1 )
!      xp1(1) = 0.
!      xp1(n1+1) = 15.         ! exp(-15)  is very small
!      
!      write(*, *) ' xp1', xp1
!      write(*, *) ' coef1', coef1
!      

      call ksampPwX(c2, p2, x02, n2, xp2(2), coef2)
      xp2(1) = 0.
      xp2(n2+1) =max( xp2(n2) * 1.5, xp1(n1+1))
!     
!      write(*, *) ' xp2', xp2
!      write(*, *) ' coef2', coef2
!     
!          merge two power functions
      call ksampPwMrg(coef1, p1, xp1, n1, 
     *                coef2, p2, xp2, n2,
     *                c,     p,   xp, nx)
!     
!      write(*,*) ' merged n c=', nx, c
!      write(*, *) ' p=', p
!      write(*,*) ' xp=', xp
!    
!    ----------------      

      do i = 1, 1000
        call ksampPw(i, c, p, xp, nx,  z, fx)
!          compute correct function value at z
        ft = cmBremF11(z) / z/(2. + 3.*up*z)**2
        call rndc(u)
        if(u .lt. ft/fx) then
           goto 200
        endif
      enddo

! -------------------------------------------------------------------
!          This rejection efficiency is as follows:
!    average trials.   standard dev.  max in 10000 gene  approx distribution
!  up = 10
!      2.52                1.94           18                 exp(-x/2.4)
!  up = 1.
!      2.805               2.23           25                 exp(-x/2.5)
!  up = 1.e-3
!      4.26                3.75           43                 exp(-x/4.0)
!  up = 1.e-6
!      almost no change from up = 1.e-3
!
! -------------------------------------------------------------------

    
!
 200  continue
       x = 3.*up*z/2.   
       x = x/(1.+x)
!      
!      write(*, *) i, sngl(x)
!      
       end
