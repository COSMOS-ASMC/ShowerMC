!
!      call kqsort(compf, idx, n)
!
!      kqsort is a quick sort routine for any data type in a 1 dimensional
!      array.  (data type may be real, integer, character, double precision
!      real or integer, half integer)
!      It can sort the data in ascending or descending order.
!      If you follow the instructions below, you will get ascending sort.
!
!      To use kqsort, you have to supply a simple integer function for each
!      array you want to sort.  The name is arbitrary, and must be like
!      
!      integer  compf
!      external compf
!      --------------------
!      integer  n
!      parameter (n = 10000)
!      real  x(n) 
!      common /abc/ x
!      ------------------- 
!      integer idx(n)
!       
!      x(*)  is computed somewhere
!      call kqsort(compf, idx, n)
!      Then, idx gets sorted order as (in default)
!      x(idx(1)) <= x(idx(2)) <= x(idx(3)).... <= x(idx(n))
!
!    ,,,
!      integer function compf(i, j)
!      integer i, j
!      --------------------
!      integer  n
!      parameter (n = 10000)
!      real  x(n) 
!      common /abc/ x
!      ------------------- 
!       if( x(i) .lt. x(j)) then
!          compf = -1     ! put  1 if you want descending sort
!       elseif( x(i) .gt.  x(j)) then
!          compf = 1      ! put -1 if you want descending sort
!       else
!          compf = 0
!       endif
!       end
!       *******************************************
!      If you want to have descending order, you may use
!      idx(n) to idx(1). However, this may lead to some confusion,
!      so you can get  descening sort directly. There are two method:
!      
!      1) After calling kqsort in a default manner, you may call
!            call ksortinv(idx, n)
!         Then, x(idx(1)) >= x(idx(2)) ... >= x(idx(n)) 
!      2) You may construct compf integer funtcion as shown by the
!         comment in the above example (reverse the sing of the function
!         value).
!
!      You may worry about the overhead of calling ksortinv,
!      but the time for it can  be negligiblly small as compared with
!      kqsort itself.
!
!      If you sort a large array (say, size >10^6)
!      many times, it  may be  better to use
!      kqsortd, kqsortr, kqsoti, kqsorth or kqsortc depending on
!      double precision real, real, integer, half integer, character
!      data. These routines don't require additonal integer function
!      like compf. 
!           call kqsort?(x, idx, n)
!      The routies are for ascending order sort; if you
!      want descending sort, use ksortinv.  
!      ********************
!
!cc       test kqsort
!      implicit none
!      integer i, j 
!      
!      external dcompf, icompf, rcompf, ccompf
!      integer dcompf, icomf, rcompf, ccompf
!      real*8  u
!
!      integer n
!      parameter (n = 1000000)
!      real*8 a( n )
!      real   b( n ) 
!      integer c( n )
!      character*9  x( 10 )
!
!      common /zzz/ a, b, c
!      common /zzzc/ x
!
!      integer idx(n)
!
!      do i = 1, n
!         call rndc(u)
!         b(i) = u
!      enddo
!      call kqsort(rcompf, idx, n)
!      do i = 1, n/2
!         j =idx(i)
!         idx(i) = idx(n-i+1)
!         idx(n-i+1) = j
!      enddo
!                                              
!      do  i=1, 10
!         write(*,*) i, b(idx(i))
!      enddo
!      do i=n-10, n
!         write(*,*) i, b(idx(i))
!      enddo
!
!
!      end
!
!
!      integer function dcompf(i, j)
!      integer i, j
!      integer n
!      parameter (n = 1000000)
!      real*8 a( n )
!      real   b( n ) 
!      integer c( n )
!      character*9  x(  10 )
!      common /zzz/ a, b, c
!      common /zzzc/ x
!
!      
!      if(a(i) .lt. a(j)) then
!         dcompf = -1
!      elseif(a(i) .eq. a(j)) then
!         dcompf = 0
!      else
!         dcompf = 1
!      endif
!      end
!      integer function rcompf(i, j)
!      integer i, j
!      integer n
!      parameter (n = 1000000)
!      real*8 a( n )
!      real   b( n ) 
!      integer c( n )
!      character*9  x(  10 )
!      common /zzz/ a, b, c
!      common /zzzc/ x
!
!      
!      if(b(i) .lt. b(j)) then
!         rcompf = -1
!      elseif(b(i) .eq. b(j)) then
!         rcompf = 0
!      else
!         rcompf = 1
!      endif
!      end
!
!      integer function icompf(i, j)
!      integer i, j
!      integer n
!      parameter (n = 1000000)
!      real*8 a( n )
!      real   b( n ) 
!      integer c( n )
!      character*9  x(  10 )
!      common /zzz/ a, b, c
!      common /zzzc/ x
!
!      
!      if(c(i) .lt. c(j)) then
!         icompf = -1
!      elseif(c(i) .eq. c(j)) then
!         icompf = 0
!      else
!         icompf = 1
!      endif
!      end
!
!      integer function ccompf(i, j)
!      integer i, j
!      integer n
!      parameter (n = 1000000)
!      real*8 a( n )
!      real   b( n ) 
!      integer c( n )
!      character*9  x( 10 )
!      common /zzz/ a, b, c
!      common /zzzc/ x
!
!      
!      if(x(i) .lt. x(j)) then
!         ccompf = -1
!      elseif(x(i) .eq. x(j)) then
!         ccompf = 0
!      else
!         ccompf = 1
!      endif
!      end
!      
      SUBROUTINE kqsort(compf, ORD,N)
!
!==============SORTS THE ARRAY A(I),I=1,2,...,N BY PUTTING THE
!   ASCENDING ORDER VECTOR IN ORD.  THAT IS ASCENDING ORDERED A
!   IS A(ORD(I)),I=1,2,...,N; DESCENDING ORDER A IS A(ORD(N-I+1)),
!   I=1,2,...,N .  THIS SORT RUNS IN TIME PROPORTIONAL TO N LOG N .
!
!
!     ACM QUICKSORT - ALGORITHM #402 - IMPLEMENTED IN FORTRAN BY
!                                 WILLIAM H. VERITY
!                                 COMPUTATION CENTER
!                                 PENNSYLVANIA STATE UNIVERSITY
!                                 UNIVERSITY PARK, PA.  16802
!
      IMPLICIT INTEGER (A-Z)
!
      DIMENSION ORD(N),POPLST(2,20)
      external compf
!c      INTEGER X,XX,Z,ZZ,Y
!c      CHARACTER *(*) A(N)
!
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
    1 ORD(I)=I
    2 IF (U1.GT.L1) GO TO 3
      RETURN
!
    3 L=L1
      U=U1
!
! PART
!
    4 P=L
      Q=U
      X=ORD(P)
      Z=ORD(Q)

!      IF (A(X).LE.A(Z)) GO TO 5
      if( compf(x, z) .le. 0 ) goto 5 
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
!
! LEFT
!
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=ORD(P)
!      IF (A(X).GE.A(XX)) GO TO 8
      if( compf(x, xx) .ge. 0) goto 8
      GO TO 6
    7 P=Q-1
      GO TO 13
!
! RIGHT
!
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=ORD(Q)
!      IF (A(Z).LE.A(ZZ)) GO TO 10
      if( compf(z, zz) .le. 0)   goto 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=ORD(P)
!
! DIST
!
!   10 IF (A(X).LE.A(Z)) GO TO 11
 10   IF ( compf(X, z) .le. 0) goto 11  
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
!   11 IF (A(X).LE.A(XX)) GO TO 12
   11 IF ( compf( x, xx) .le. 0) goto 12
      XX=X
      IX=P
!   12 IF (A(Z).GE.A(ZZ)) GO TO 6
   12 IF ( compf(z, zz) .ge. 0) goto 6
      ZZ=Z
      IZ=Q
      GO TO 6
!
! OUT
!
   13 CONTINUE
!      IF (.NOT.(P.NE.IX.AND.A(X).NE.A(XX))) GO TO 14
      IF (.NOT.(P.NE.IX.AND. compf(X, xx) .ne. 0) ) goto 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
!      IF (.NOT.(Q.NE.IZ.AND.A(Z).NE.A(ZZ))) GO TO 15
      IF (.NOT.(Q.NE.IZ.AND. compf(z, zz) .ne. 0) ) goto 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18
!
! START RECURSIVE CALL
!
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
!
! POP BACK UP IN THE RECURSION LIST
!
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
!
! END SORT
! END QSORT
!
      END
