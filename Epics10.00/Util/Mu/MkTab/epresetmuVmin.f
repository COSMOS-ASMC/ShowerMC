      subroutine epresetmuVmin(vmin)
      implicit none
!        read vmin from /tmp/muVmin  and
!        reset vmin=Et/Emu for muon brems, pair or nuc. int.
!
      real*8 vmin        ! output. vmin= Et/Emu

      integer icon, io

      io = 11

      call copenf(io, './muVmin', icon)
      if(icon .ne. 0) then
         call cerrorMsg('./muVmin not exit', 0)
      endif
      read(io, *) vmin
      close(io)
      end

