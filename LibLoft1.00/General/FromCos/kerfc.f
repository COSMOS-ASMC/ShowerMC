!           test kerfc
!      implicit none
!      real*8 x, kerfc, kerf
!      do x= 0., 20., 0.5
!         write(*, *) x, " ", kerfc(x), ' ', 1.d0-kerf(x)
!      enddo
!      end
!     ***********************
      real*8 function kerfc(x)
!     *********************** 
! 
!     complementary error function 
!     = 1 - kerf(x).
!     Note. for x > 4, 1 - kerf(x) gives
!      poor result. Therefore we must use
!      kerfc directly.
!
!    Adapted from Mori's code in Maruzen book.
!
      implicit none

      real*8 x
      integer nm, nx, na
      parameter (nm = 5, nx = 13, na = 5)

      real*8 cm(0:nm),cx(nx),cq(nx),ca(0:na)
      real*8 sqrtpi, invpi, norm, cxi
      parameter (sqrtpi = 1.772453850905516d0)
      parameter (norm = 2.d0 / sqrtpi)

      parameter (invpi = 1. / 3.141592653589793d0)
      parameter (cxi = 4.d0 / invpi)
      real*8 xv, y, v
      integer i

      data cm /
     *   0.1000000000000000d+01,
     *  -0.3333333333333333d+00,
     *   0.1000000000000000d+00,
     *  -0.2380952380952381d-01,
     *   0.4629629629629630d-02,
     *  -0.7575757575757575d-03 /
!
      data cx /
     *   0.7788007830714048d+00,
     *   0.3678794411714423d+00,
     *   0.1053992245618643d+00,
     *   0.1831563888873418d-01,
     *   0.1930454136227709d-02,
     *   0.1234098040866796d-03,
     *   0.4785117392129009d-05,
     *   0.1125351747192591d-06,
     *   0.1605228055185612d-08,
     *   0.1388794386496402d-10,
     *   0.7287724095819692d-13,
     *   0.2319522830243569d-15,
     *   0.4477732441718302d-18 /
      data cq /
     *   0.2500000000000000d+00,
     *   0.1000000000000000d+01,
     *   0.2250000000000000d+01,
     *   0.4000000000000000d+01,
     *   0.6250000000000000d+01,
     *   0.9000000000000000d+01,
     *   0.1225000000000000d+02,
     *   0.1600000000000000d+02,
     *   0.2025000000000000d+02,
     *   0.2500000000000000d+02,
     *   0.3025000000000000d+02,
     *   0.3600000000000000d+02,
     *   0.4225000000000000d+02 /
!
      data ca /
     *   0.1000000000000000d+01,
     *  -0.1000000000000000d+01,
     *   0.3000000000000000d+01,
     *  -0.1500000000000000d+02,
     *   0.1050000000000000d+03,
     *  -0.9449999999999999d+03 /
!

      xv = abs(x)

      if (xv .le. 0.1d0) then
         y = xv**2
         v = cm(nm)
         do i = nm - 1, 0, -1
            v = cm(i) + y * v
         enddo
         kerfc = 1.d0  - norm * xv * v
      elseif(xv .le. 100.0d0) then
        y = xv**2
        v = 1 / (2 * y)
        do i = 1, nx
           v = v + cx(i) / (cq(i) + y)
        enddo
        v = invpi * xv * exp(-y) * v
        if (xv .lt. 6.0d0) then
           v = v - 2 / (exp(cxi * xv) - 1)
        end if
        kerfc = v
      else
        y = 2 * xv**2
        v = ca(na)
        do i = na - 1, 0, -1
          v = ca(i) + y * v
        enddo
        v = exp(-xv**2) / (sqrtpi * xv) * v
        kerfc = v
      endif

      if (x .lt. 0.) then
        kerfc = 2.d0  - kerfc
      endif
      end
