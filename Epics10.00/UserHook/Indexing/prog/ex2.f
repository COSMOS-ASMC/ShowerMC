!      This is for config (CALET). so put config in FirstInput
      module modIndexing
!           SciFi for X,Y 
      character(len=*),
     *  parameter:: SciFiSpec="|sfx1sheet,sfy1sheet|  SciFi"
      integer:: SciFiShape(3)   !  if | | is used, rank must be 2+1 
      integer,allocatable::  SciFiIdx(:,:,:)
      integer,pointer::SciFinum(:,:)
!         reshaped index; Fortran oriented version
      integer:: FSciFiShape(3)
      integer,allocatable::  FSciFiIdx(:,:,:)
      end       module modIndexing

      subroutine epSetIdx
      use modIndexing
      use modGetIndex
      implicit none
      integer:: i,j,  n
      integer:: Nlayer, NSciFi

      call epCountSubdTree(SciFiSpec, SciFiShape, n, SciFinum)
      write(0,*) 'SciFinum =',SciFinum
      Nlayer =  ScifiShape(2)
      NSciFi =  ScifiShape(3)
                                   ! =2
      write(0,*) ' ScifiShape=', SciFiShape(1),  Nlayer, NSciFi
      
      allocate( SciFiIdx(2, Nlayer, NSciFi) )

      call epGetIndex(SciFiSpec, SciFiShape, SciFiIdx, SciFinum)

      do i = 1, Nlayer
         do j = 1, NSciFi
            write(0,'(a, i4, a, i4, a, i6)')
     *      ' ScifiX # ',j, ' in layer ',i, ' has  comp#=',
     *      SciFiIdx(1, i,j)
         enddo
      enddo

      do i = 1, Nlayer
         do j = 1, NSciFi
            write(0,'(a, i4, a, i4, a, i6)')
     *      ' ScifiY # ',j, ' in layer ',i, ' has  comp#=',
     *      SciFiIdx(2, i,j)
         enddo
      enddo
!           reshape; Fortran oriented array
      forall(i=1:3) FSciFiShape(4-i) =SciFiShape(i)

      allocate( FSciFiIdx(NSciFi, Nlayer, 2) )
      FSciFiIdx =
     *      reshape(SciFiIdx, FSciFiShape, order=(/3,2,1/))
      do j = 1, Nlayer
         do i = 1, NSciFi
            write(0,'(a, i4, a, i4, a, i6)')
     *      'F ScifiX # ',i, ' in layer ',j, ' has  comp#=',
     *           FSciFiIdx(i,j,1)
         enddo
      enddo
!           Y
      do j = 1, Nlayer
         do i = 1, NSciFi
            write(0,'(a, i4, a, i4, a, i6)')
     *      'F ScifiY # ',i, ' in layer ',j, ' has  comp#=',
     *           FSciFiIdx(i,j,2)
         enddo
      enddo

      end subroutine epSetIdx
