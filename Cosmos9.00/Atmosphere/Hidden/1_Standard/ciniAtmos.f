      subroutine ciniAtmos
      implicit none
      integer::atmosmodel
!      alloate tables      
      call callocAtmosTbl
!          read segmented atmosphere data
      call creadAtmosD
!          manipulate data
      call catmosCnst1
      call cqAtmosModel(atmosmodel)
      write(0,*) ' Atmosphere model # is ', atmosmodel
      end   subroutine ciniAtmos
