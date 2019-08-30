#undef  ABSOFT
#ifdef  NEXT486
#define  ABSOFT
#elif  defined PCLINUX
#define  ABSOFT
#elif  defined MACOSX
#define  ABSOFT
#endif
!
!           manipulate bit
!
      logical function kbitest( i, bit )
      implicit none
      integer i  ! input.  it's bit-th bit is tested.
      integer bit ! input.  (1<- bit <=31) 
!       function value. .true. if bit-th bit of i is on.
!                       else .false.
!                       if bit is outside of the range,
!                       .false. is returned.
!
      integer j, k
#ifdef  ABSOFT
      kbitest = btest(i, bit-1)
#else
      
      if( bit .le. 0 .or. bit .ge. 32) then
         kbitest = .false.
      else
         j = lshift(i, 32-bit)
         k = rshift(j, 31)
         kbitest = k .ne. 0
      endif
#endif
      end
      subroutine ksetbit( i, bit )
      implicit none
      integer i  !  input/output.
      integer bit  !  input.  bit-th position of i is made to be on.
!
      logical kbitest
      integer j, k
#ifdef ABSOFT
      i = ibset(i, bit-1)
#else
      if(bit .gt. 0 .and. bit .lt. 32) then
         if( .not. kbitest(i, bit) ) then
            i = i + 2**(bit-1)
         endif
      endif
#endif
      end
      subroutine krsetbit( i, bit )
      implicit none
      integer i  !  input/output.
      integer bit  !  input.  bit-th position of i is made to be off.
!
      logical kbitest
      integer j, k
#ifdef ABSOFT
      i = ibclr(i, bit-1)
#else
      if(bit .gt. 0 .and. bit .lt. 32) then
         if( kbitest(i, bit) ) then
            i = i - 2**(bit-1)
         endif
      endif
#endif
      end
      subroutine krevbit( i, bit )
!       reverse bit-th bit of i.
      implicit none
      integer i  !  input/output.
      integer bit  !  input.  bit-th position of i is made to be reverted
!
      logical kbitest
      integer j, k
#ifdef ABSOFT
      if(kbitest(i, bit)) then
         i = ibclr(i, bit-1)
      else
         i = ibset(i, bit-1)
      endif
#else
      if(bit .gt. 0 .and. bit .lt. 32) then
         if( kbitest(i, bit) ) then
            i = i - 2**(bit-1)
         else
            i = i + 2**(bit-1)
         endif
      endif
#endif
      end
      
!       test kbitest.
!
!      integer i, j
!      logical kbitest
!
!      do while ( .true. )
!         read(*, *) i, j
!         write(*,*) i, j, kbitest(i, j)
!         call ksetbit(i, j)
!         write(*, *) i, kbitest(i, j)
!         call krsetbit(i, j)
!         write(*, *) i, kbitest(i, j)
!         call krsetbit(i, j)         
!         write(*, *) i, kbitest(i, j)
!      enddo
!      end




