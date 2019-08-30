!          make random number seed using timer, hostname, and process number.
      subroutine cmkSeed(dummy, seed)
      implicit none

      integer dummy    ! input. some integer used to make seed
      integer seed(2)  ! output.  two integers seed to be used as
                       !          call rndir(seed)




      integer leng, i, j, kgetpid, kgettime
      integer k
      character*24 hostn



      seed(1)= kgettime(dummy) + dummy
      seed(2)= kgetpid(dummy)      ! this is really dummy 
!
!         add a max of 4 last char. of hostname after converting to int.
!
      call cgetHost(leng, hostn) ! get hostname
      i = max(1, leng-3)

      do  j = i, leng
         seed(2) = seed(2) + ichar( hostn(j:j) )
         seed(1) = seed(1) + ichar( hostn(j:j) )
      enddo
!
!        next is from rnd1u: to avoid that the 
!        first random number is biased to large or
!        small number systematically
!
      k = seed(1)/53668
      seed(1) = 40014*(seed(1) - k*53668) - k*12211
      if(seed(1) .lt. 0) seed(1)=seed(1)+2147483563
      k = seed(2)/52774
      seed(2) = 40692*(seed(2) - k*52774) - k* 3791
      if(seed(2) .lt. 0) seed(2)= seed(2)+2147483399

      end
