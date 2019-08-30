!      implicit none
!      real*8 kzeta2a, a
!      do a=0., 1., 0.02
!         write(*, *) sngl(a), sngl(kzeta2a(a))
!      enddo
!      end
!     ****************************************************************
!     *                                                              *
!     * kzeta2a: compute sum of a**m/m**2 from m=1 to inf ( 0<=a<=1) *
!     
!            a       a2      a3       a4
!           ---  +  ---  +  ----  +  ---
!              2       2       2        2
!            1        2       3        4
!     ****************************************************************
!
!    usage: real*8 kzeta2a     
!            f = kzeta2a(a)
!
!  method:   use the series as it is if a < .5 to get sum with
!            relative error 1/1000 else
!            use polinomial approximation to result obtained by
!            equivalent integral ( 0 to -ln(1-a) of y/(exp(y) -1) )
!
!            zeta2a(1.) = zeta(2) =pi**2/6 = 1.644934
!
! note: accuracy is only a few to several digits.
!
      real*8  function kzeta2a(a)
      implicit none
      real*8 a
!
      real*8 c(7)/
     * 29.55194, -260.5342, 955.5774, -1843.917, 1984.555, -1129.490,
     * 265.9021/

      real*8 s, r, ak, tmp
      integer k, j
!
      s=0.
      if(a .lt. .5) then
          r=1.
          if(a .eq. 0.) r=0.
          ak=a
          k=1
          do   while (r .gt. 1.e-3 )
              tmp= ak/k**2
              s=s+tmp
              r=abs(tmp/s)
              ak=ak*a
              k=k+1
          enddo
      else
           do   j=7, 1, -1
             s= s *a + c(j)
           enddo
      endif
      kzeta2a=s
      end

