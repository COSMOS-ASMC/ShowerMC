!c         test copenf and cskipComment
!     implicit none
!     character*80 string, fn

!     data fn/'../Data/Primary/sample.dat'/

!     integer icon
!     call copenf(TempDev, fn, icon)
!     write(*, *) " icon=", icon
!     call cskipComment(TempDev, icon)
!     write(*, *) ' icon= ', icon
!     read(TempDev, '(a)') string
!     write(*, *) string
!     end
!     kskip_commnet *****************************************************
!     skip data until a line starting with "#------------   "
!
       subroutine cskipComment(io, icon)
!         io: integer. input.logical file #
!       icon: integer. output. 0 --> ok
!                              1 --> no such line or error
!
       implicit none
       integer io, icon

       character*16 headline
       integer ios
       logical ok

       ok = .false.
       icon = 0
       do while(.not. ok)
          read(io, '(a)', iostat=ios) headline
          if(ios .ne. 0) then
              icon = 1
              ok = .true.
              write(*, *) ' no separator'
          else
              ok =headline(1:16) .eq. '#---------------'
          endif
       enddo
       end
