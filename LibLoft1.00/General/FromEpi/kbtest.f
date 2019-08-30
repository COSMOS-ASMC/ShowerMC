!      implicit none
!      integer i, j, kbclr, kbset
!      logical kbtest
!
!      do while( .true. )
!         i = 12345
!         read(*,*) i,j
!         i = kbset(i, j)
!         write(*,*) ' set ',j,'-th bit  ', i
!         write(*,*) ' test it=', kbtest(i,j)
!         i = kbclr(i, j)
!         write(*,*) ' reset ',j,'-th bit  ', i
!         write(*,*) ' test it=', kbtest(i,j)
!      enddo
!      end

      integer function kbset(bitpat, bitpos)
      implicit none
      integer bitpat  ! input. see kbtest
      integer bitpos  ! input.  note lsb  is bit  pos 1.

      if(bitpos .gt. 32 .or. bitpos .lt. 1) then
         kbset = bitpat
      else
         kbset = ibset(bitpat, bitpos-1)
      endif
      end
!
      integer function kbclr(bitpat, bitpos)
      implicit none
      integer bitpat   ! input.  see  kbtest
      integer bitpos   ! input.  

      if(bitpos .gt. 32 .or.  bitpos .lt. 1) then
         kbclr = bitpat
      else
         kbclr = ibclr(bitpat, bitpos-1)
      endif
      end

      logical  function kbtest(bitpat, bitpos)
      implicit none
      integer bitpat  ! bitpattern to be examined
      integer bitpos  ! examine if bitpos-th bit of bitpat
                      ! is on or not.  if  on,  function value
                      ! becomes .true. else .false.
                      ! the LSB  is bit  position 1.
                      ! bitpos must be 1-32.  For other values,  function
                      ! value becomes .false.

      if( bitpos .gt. 32 .or.  bitpos .lt. 1) then
         kbtest = .false.
      else
         kbtest = btest(bitpat, bitpos-1)
      endif
      end

           
