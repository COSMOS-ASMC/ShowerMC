      subroutine  fort
      real*8   xyz
      integer  yyy
      common /abc/ xyz, yyy

      integer   iab
      common /ppp/ iab

      write(*,*) xyz, yyy, iab
      end

