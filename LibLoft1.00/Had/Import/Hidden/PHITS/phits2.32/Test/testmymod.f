      program main
      use mymodx
      implicit none
      call sub1
      call sub3
      call sub2
      end program main

      subroutine sub2
      write(0,*) ' this is sub2 not in module'
      end
