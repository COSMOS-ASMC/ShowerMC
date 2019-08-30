!       fix  binorasc:  1 for ascii output (myhist1.hist)
!                       2 for binary //    (myhist1.chist)
!        use   call kwhists0(0)  for integral from -inf.
!           or call kwhists0(1)  fro integral from +inf.
  !
program main
  use modHistogram
  use modHistogram1
  implicit none

  integer binorasc/1/
  integer nav, nsig, npw
  parameter (nav=2, nsig=5, npw=3 )

  type(histogram1)  h(nav, nsig)
  type(histogram1)  k(npw)
  real*8 av, sig, pw
  save h
  real*8 x
  integer i, j, m, fno
  integer nc
  character*48 dirstr
  character*38 key
  real*8 xxx(500), yyy(500)
!  external kwhistIxy
  !  integer  nnn, kwhistIxy
  integer  nnn

  fno= 3
  if(binorasc .eq. 1) then
     open(fno, file='mytest1.hist', form='formatted')
  else
     open(fno, file='mytest1.chist', form='unformatted')
  endif

  call kwhistso(binorasc)  ! bin/asc write.  this is common to all.

!       minimum calls
  pw = 0.8
  do i = 1, npw
     pw = pw + 0.2
!              init.
     call kwhisti(k(i), 1.5, 0.1, 30, b'01111' )
!              clear
     call kwhistc(k(i))
     do j= 1, 1000000
        call rndc(x)
        x = x**(-pw)
        !               take histo
        call kwhist( k(i), sngl(x), 1.0 )
     enddo
  enddo
!             make arg 0 for integral from -inf. 
!                      1 for //            +inf.
  call kwhists0(0)
  !        output 
  do i =1, npw
     call kwhists(k(i), 0.)
     call kwhistp(k(i), fno )
  enddo

! ++++++++++++++++++++++++++++++
!    some standard
!
  av=0.
  do i = 1, nav
     av = av + 2. 
     sig = 0.2 
     do j = 1, nsig
        sig = sig + 0.2
        !             init.
        call kwhisti(h(i, j), sngl(av-7*sig), sngl(av+7*sig),   &
             500,  b'10000')
        !             clear
        call kwhistc(h(i, j))
        do m = 1, 1000000
           call kgauss(av, sig, x)
           !                 take histo
           call kwhist(h(i,j), sngl(x), 1.0 )
        enddo
!             give additional info.
        call kwhistai(h(i,j),  &
                "Test Gaussian dist.",  &
                "gauss", "event", .false., 0.,  &
                "x", "m") 
!             make key for diff. parameers
        write(key,'(1p2E11.3)' ) av, sig
!              inform it
        call kwhistid(h(i,j), key)
!              make directory: maindir/gauss/{av1,av2}
        write(dirstr,'("av",i2,"/")')  i
!                   next two is better to shrink the string length
!                   but only 3rd line can be  ok. (white blank will 
!                   be eliminated inside).
!            call kseblk(dirstr,"|", nc)
!            call kwhistdir(h(i,j), dirstr(1:nc))
        call kwhistdir(h(i,j), dirstr)
     enddo
  enddo
!         output.
  do i = 1, nav
     do j = 1, nsig
        call kwhists(h(i,j), 0.)
        call kwhistp(h(i,j), fno)
     enddo
  enddo
!            get (x,y) of integral dist.
  nnn=kwhistIxy(h(2,2), xxx, yyy, 500)
  write(0,*) ' nnn-', nnn
  do i = 1, nnn
     write(0,*) xxx(i), yyy(i)
  enddo
end program main
