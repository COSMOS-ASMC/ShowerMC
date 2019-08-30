!      This is for config (CALET). so put config in FirstInput
      module modIndexing
!           SciFi for X
      character(len=*),parameter:: SciFiXSpec="sfx1sheet sfx"
      integer:: SciFiXShape(2)
      integer,allocatable::  SciFiXIdx(:,:)
      integer,pointer::SciFiXnum(:,:)
!         reshaped index; Fortran oriented version
      integer:: FSciFiXShape(2)
      integer,allocatable::  FSciFiXIdx(:,:)
!!!!!!!!!!!!!
!          SciFi for Y 
      character(len=*),parameter:: SciFiYSpec="sfy1sheet sfy"
      integer:: SciFiYShape(2)
      integer,allocatable::  SciFiYIdx(:,:)
      integer,pointer::SciFiYnum(:,:)
!!!!!!!!!!!!!
!       For  X,Y combined idex      
      integer:: SciFiShape(3)
      integer,allocatable::  SciFiIdx(:,:,:)
!
!          other SciFi
      character(len=*),parameter:: specOptX ="optx14u SciFi"
      integer:: OptXShape(2)
      integer,allocatable:: OptXIdx(:,:)
      integer,pointer:: OptXnum(:,:)

      character(len=*),parameter:: specOptY ="opty14u SciFi"
      integer:: OptYShape(2) 
      integer,allocatable:: OptYIdx(:,:)
      integer,pointer:: OptYnum(:,:)
      end       module modIndexing

      subroutine epSetIdx
      use modIndexing
      use modGetIndex
      implicit none
      integer:: i,j,  n
      integer:: Nlayer, NSciFi

      call epCountSubdTree(SciFiXSpec, SciFiXShape, n, SciFiXnum,
     *  target="subd"  )
      write(0,*) 'SciFiXnum =',SciFiXnum
      Nlayer =  ScifiXShape(1)
      NSciFi =  ScifiXShape(2)
      write(0,*) ' ScifiXShape=', Nlayer, NSciFi
      
      allocate( SciFiXIdx(Nlayer, NSciFi) )

      call epGetIndex(SciFiXSpec, SciFiXShape, SciFiXIdx, SciFiXnum,
     *  target="subd")

      do i = 1, Nlayer
         do j = 1, NSciFi
            write(0,'(a, i4, a, i4, a, i6)')
     *      ' ScifiX # ',j, ' in layer ',i, ' has  comp#=',
     *      SciFiXIdx(i,j)
         enddo
      enddo
!           reshape; Fortran oriented array
      forall(i=1:2) FSciFiXShape(3-i) =SciFiXShape(i)

      allocate( FSciFiXIdx(NSciFi, Nlayer) )
      FSciFiXIdx =
     *      reshape(SciFiXIdx, FSciFiXShape, order=(/2,1/))
      do j = 1, Nlayer
         do i = 1, NSciFi
            write(0,'(a, i4, a, i4, a, i6)')
     *      'F ScifiX # ',i, ' in layer ',j, ' has  comp#=',
     *           FSciFiXIdx(i,j)
         enddo
      enddo
!           Y
      call epCountSubdTree(SciFiYSpec, SciFiYShape, n, SciFiYnum,
     *   target="subd")
      write(0,*) 'SciFiYnum=',SciFiYnum

      Nlayer =  ScifiYShape(1)
      NSciFi =  ScifiYShape(2)
      write(0,*) ' ScifiYShape=', Nlayer, NSciFi
      
      allocate( SciFiYIdx(Nlayer, NSciFi) )
      call epGetIndex(SciFiYSpec, SciFiYShape, SciFiYIdx, SciFiYnum,
     *   target="subd")
!         combine two 
      SciFiShape(1) = 2
      SciFiShape(2:3) = SciFiXShape(1:2)

      allocate( SciFiIdx(2, Nlayer, NScifi) )
      SciFiIdx(1,:,:) = SciFiXIdx(:,:)
      SciFiIdx(2,:,:) = SciFiYIdx(:,:)

      deallocate( SciFiXIdx )
      deallocate( SciFiYIdx )

      write(0,*) ' ' 
      do i = 1, Nlayer
         do j = 1, NSciFi
            write(0,'(a, i4, a, i4, a, i6)')
     *           ' ScifiX # ',j, ' in layer ',i, ' has  comp#=',
     *           SciFiIdx(1, i,j)
         enddo
      enddo
      write(0,*) ' ' 
      do i = 1, Nlayer
         do j = 1, NSciFi
            write(0,'(a, i4, a, i4, a, i6)')
     *      ' ScifiY # ',j, ' in layer ',i, ' has  comp#=',
     *           SciFiIdx(2, i,j)
         enddo
      enddo
      end subroutine epSetIdx
