      subroutine ciniAtmos
      implicit none
      integer:: atmosmodel
!          read segmented atmosphere data
      call creadAtmosD
!          manipulate data
      call catmosCnst1
      call catmosCnst2
      call cqAtmosModel(atmosmodel)
      write(0,*) ' Atmosphere model # is ', atmosmodel
      end   subroutine ciniAtmos
      
