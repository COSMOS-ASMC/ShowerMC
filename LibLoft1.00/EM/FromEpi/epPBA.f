      subroutine epBremAng(e, m, eg, z, teta)
      implicit none
      real*8  e ! input. energy of electron/positron GeV
      real*8  m ! input. electron  mass  GeV.
      real*8 eg ! input. brems gamma energy.  GeV.
      real*8  z ! input. matter z.
      real*8  teta ! output. sampled angle of photons relative to
                   !       the inciedent electron
      real*8 d, pi, hpi, maxu, u
      real*8 a/0.625d0/

      
!
      parameter( pi=3.14159265d0, hpi= pi/2)

      d = 0.13d0*(0.8d0 + 1.3d0/z)*
     *         (100.d0+ 1.d0/e) * (1.+ eg/e)

!
!      maxu = e*hpi/m
      maxu = e*pi/m   ! 2018,Apr.3
      do while (.true.)
         call epPBAng(a, d, u)
         if(u .lt. maxu) goto 10
      enddo
 10   continue
!
      teta = u*m/e
      end

      subroutine epPairAng(e, m,  teta)
      implicit none
!          samples polar angle of an electron or positron 
!        from  pair creation.  It must be a smaller eneregy
!        one.
!        The angular distribution is approximated as
!      {  uexp(-au) + duexp(-3au) } du
!        where a = 0.625 and d=27.0. This is  employed in
!        G3
      real*8 e ! input. energy of e+/e_ (smaller energy)
      real*8 m ! input. mass of the electron. must be in the same unit as
               !         e

      real*8 teta ! output.  sampled angle in radian. < pi/2
!                          
      real*8 a/0.625d0/, d/27.0d0/, pi, hpi, maxu, u
      parameter( pi=3.14159265d0, hpi= pi/2)
!
      maxu = e*hpi/m
      do while (.true.)
         call epPBAng(a, d, u)
         if(u .lt. maxu) goto 10
      enddo
 10   continue
!
      teta = u*m/e
      end
      subroutine epPBAng(a, d, u)
      implicit none
!      sample u from
!       {  uexp(-au) + duexp(-3au) } du
!
      real*8 a  ! input. > 0. see above.
      real*8 d  ! iput.  >= 0  see above.
!
      real*8 u  ! output. sampled variable u.
!         the integration gives
!
!             1/a^2 + d/9/a^2 so that the weight is
!          1 : d/9
!        9: d =   9/(9+d):  d/(9+d)
!
      real*8  x, x1, x2
      call rndc(x1)
      call rndc(x2)
      u = -log(x1*x2) ! this obeys uexp(-u)du
      call rndc(x)
      if(x .lt. 9./(9.+d)) then
!          use  u exp(-au)du
         u = u/a
      else
!          use u exp(-3au) du
         u = u/3.d0/a
      endif
      end
      


