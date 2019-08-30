!      This is for config (CALET). so put config in FirstInput
!  to make this work, you have to make
!     optx14u and opty14u as a world.  NOT implicit world
!     but add _w to the last "horse" like "horse_w"

      module modIndexing
!          other SciFi
      character(len=*),parameter:: OptXspec ="optx14u SciFi"
      integer:: OptXShape(2)
      integer,allocatable:: OptXIdx(:,:)
      integer,pointer:: OptXnum(:,:)

      character(len=*),parameter:: OptYspec ="opty14u SciFi"
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

      call epCountSubdTree(OptXspec,  OptXShape, n, OptXnum)
      write(0,*) 'OptXnum =',OptXnum
      Nlayer =  OptXShape(1)
      NSciFi =  OptXShape(2)
      write(0,*) ' OptX Shape=', Nlayer, NSciFi
      
      allocate( OptXIdx(Nlayer, NSciFi) )

      call epGetIndex(OptXspec, OptXShape, OptXIdx, OptXnum)

      do i = 1, Nlayer
         do j = 1, NSciFi
            write(0,'(a, i4, a, i4, a, i6)')
     *      ' ScifiX # ',j, ' in layer ',i, ' has  comp#=',
     *      OptXIdx(i,j)
         enddo
      enddo
!           Y
      call epCountSubdTree(OptYspec, OptYShape, n, OptYnum)
      write(0,*) 'OptYnum=',OptYnum(:,:)

      Nlayer =  OptYShape(1)
      NSciFi =  OptYShape(2)
      write(0,*) ' OptYshape=', Nlayer, NSciFi
      
      allocate( OptYIdx(Nlayer, NSciFi) )
      call epGetIndex(OptYSpec, OptYShape, OptYIdx, OptYnum)

      do i = 1, Nlayer
         do j = 1, NSciFi
            write(0,'(a, i4, a, i4, a, i6)')
     *      ' ScifiY # ',j, ' in layer ',i, ' has  comp#=',
     *      OptYIdx(i,j)
         enddo
      enddo

      end subroutine epSetIdx
