!  Give .hyb datfile.  It may or may not  contain  a blank line as an event separator.
! complile: e.g.  ifort integratedEdx.f -L$COSMOSTOP/lib/MacIFC -lcosmos
! usage echo  xxx.hyb | ./a.out
      implicit none
      integer np, nc, i
      parameter(np=100, nc=100)
      real*8 x(np), y(np)
      real*8 dmy
      integer npa, idmy, eventno
      character*130 hybfile, input
      hybfile = ' '
      read(*, '(a)')  hybfile
      open(11, file=hybfile)
      npa = 0
      eventno= 0
      do while(.true.)
         read(11, '(a)', end=999) input
         if( input(1:1) .eq. 'h' ) then
            if(npa .gt. 0) then
               eventno = eventno + 1
               call integrate(eventno, x, y, npa)
               npa = 0
            endif
            continue
         elseif( input  .eq.  ' ' ) then
            continue
         else
            npa = npa + 1
!            read(input,
!     *           '(1x, i4, f7.1, f7.1, 2f6.3,5e11.3,e11.3)')
!     *        idmy, x(npa), dmy, dmy, dmy, dmy, dmy, dmy, dmy,
!     *        dmy, y(npa)
            read(input(2:130),*) 
     *        idmy, x(npa), dmy, dmy, dmy, dmy, dmy, dmy, dmy,
     *         dmy, y(npa)
         endif
      enddo   
 999  continue
      if(npa .gt. 0) then
         eventno = eventno + 1
         call integrate(eventno, x, y, npa)
      endif
      close(11)
      end
      subroutine integrate(eventno, x, y, npa)
      implicit none
      integer npa, eventno
      integer np, nc
      parameter(np=100, nc=100)
      real*8 x(np), y(np)
      real*8 ans, coef(nc, 3)
      real*8 a, b, s
       integer i


      a = x(1)
      b = x(npa)
      call ktrpzIntT2(y, 1, npa, x, 1, a, b, ans)
      call kcsplCoef(x, y, npa, coef, nc)
      call kcsplInteg(x, y, npa,  coef, nc, a, b, s)
      write(*,'(i5, 1p 2g14.3, i4)') eventno, ans, s, npa
      end

