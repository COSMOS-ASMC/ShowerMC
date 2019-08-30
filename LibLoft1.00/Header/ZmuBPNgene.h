      common /muBPNgene/ media, Emu, normby
       type(epmedia)::   media
      real*8 Emu
      integer normby
!
!      normby: 1-->  mb/molecule normalization
!              2-->  /X0  normalization
!              This must be set before calling
!              muon brems/pair/nuclear  interaction
!              routines 
