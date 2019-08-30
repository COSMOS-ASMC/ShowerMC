!      This is for config (CALET). so put config in FirstInput
      module modIndexing
!           SciFi for X,Y 
      character(len=*),
     *  parameter:: SciFiSpec="(sfx1sheet,sfy1sheet)  SciFi"
      integer:: SciFiShape(2)   !  if | | is used, rank must be 2+1 
      integer,allocatable::  SciFiIdx(:,:)
      integer,pointer::SciFinum(:,:)
!         reshaped index; Fortran oriented version
      integer:: FSciFiShape(2)
      integer,allocatable::  FSciFiIdx(:,:)

      end       module modIndexing

      subroutine epSetIdx
      use modIndexing
      use modGetIndex
      implicit none
      integer:: i,j,  n
      integer:: Nlayer, NSciFi

      call epCountSubdTree(SciFiSpec, SciFiShape, n, SciFinum)
      write(0,*) 'SciFinum =',SciFinum
      Nlayer =  ScifiShape(1)
      NSciFi =  ScifiShape(2)
                                
      write(0,*) ' ScifiShape=',  Nlayer, NSciFi
      
      allocate( SciFiIdx( Nlayer, NSciFi) )

      call epGetIndex(SciFiSpec, SciFiShape, SciFiIdx, SciFinum)

      do i = 1, Nlayer
         do j = 1, SciFinum(1,1)
            write(0,'(a, i4, a, i4, a, i6)')
     *      ' Scifi(X) # ',j, ' in layer ',i, ' has  comp#=',
     *      SciFiIdx(i,j)
         enddo
      enddo

      do i = 1, Nlayer
         do j = 1, SciFinum(1,2)
            write(0,'(a, i4, a, i4, a, i6)')
     *      ' Scifi(Y) # ',j, ' in layer ',i, ' has  comp#=',
     *      SciFiIdx(i,j+SciFinum(1,1))
         enddo
      enddo
!           reshape; Fortran oriented array
      forall(i=1:2) FSciFiShape(3-i) =SciFiShape(i)

      allocate( FSciFiIdx(NSciFi, Nlayer) )
      FSciFiIdx =
     *      reshape(SciFiIdx, FSciFiShape, order=(/2,1/))
      do j = 1, Nlayer
         do i = 1, SciFinum(1,1)
            write(0,'(a, i4, a, i4, a, i6)')
     *      'F Scifi(X) # ',i, ' in layer ',j, ' has  comp#=',
     *           FSciFiIdx(i,j)
         enddo
      enddo
!           Y
      do j = 1, Nlayer
         do i = 1, SciFinum(1,2)
            write(0,'(a, i4, a, i4, a, i6)')
     *      'F Scifi(Y) # ',i, ' in layer ',j, ' has  comp#=',
     *           FSciFiIdx(i+SciFinum(1,1),j)
         enddo
      enddo

      end subroutine epSetIdx
