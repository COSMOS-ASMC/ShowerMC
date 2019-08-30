!
!          Trapezoidal integral of a given table data.
!          of which x interval is not unique.
      subroutine ktrpzIntT2(t, intv, n, xt, intvx, a, b, ans)
      implicit none
      integer intv  ! input. see below
      integer n     ! input. number of data values
      real*8 t(intv, n)  ! input.  t(1, 1), t(1, 2), .. t(1, n) are used.
                         ! function values at xt(1,1), xt(1,2), ..
      integer intvx     ! see below
      real*8 xt(intvx, n)  ! input.   must be xt(1,i) < xt(1,i+1)
      real*8 a   ! input. lower value of the integral region.
      real*8 b   ! input. upper value of the integral region.  a<= b.
!           Note:  table values outside of the given xt are assumed to be 0.
      real*8 ans ! output.  integral value

      integer i, i1, i2
      real*8 aa, bb, fa, fb
      


      ans = 0.
      aa = max(a, xt(1, 1))
      bb = min(b, xt(1, n))  
      do i =1, n
         if(aa .le. xt(1,i)) goto 10
      enddo
!         never come here

 10   continue
      i1 = i

      do i = n, 1, -1
         if(bb .ge. xt(1,i)) goto 20
      enddo
!        never come here

 20   continue
      i2 = i
!     
      if(i1 .ne. 1) then
         fa =
     *        (t(1,i1)- t(1,i1-1))/ (xt(1,i1) - xt(1,i1-1)) *
     *        (aa - xt(1,i1-1))  + t(1,i1-1)
         ans = ans + (xt(1,i1)- aa) * (fa + t(1,i1))/2
      endif
      if(i2 .ne. n) then
         fb = (t(1,i2+1) - t(1, i2))/ (xt(1,i2+1) - xt(1,i2)) *
     *        (bb -xt(1,i2)) + t(1,i2)          
         ans = ans+ (bb - xt(1,i2)) * (t(1,i2) + fb)/2
      endif
      do i = i1, i2-1
         ans = ans + (xt(1,i+1) - xt(1, i)) * (t(1,i) + t(1,i+1))/2
      enddo
      end
