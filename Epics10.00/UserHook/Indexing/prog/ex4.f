!       For testing odd detector attachment 
!          configB is used;  to be specified in FirstInput
      module modIndexing
      character(len=*),parameter:: 
     *  Spec="(sheet,oddx4) CHDx"   ! only 1 oddx4 must  accompany  sheet
!          next  can be used without changing other part
!     *  Spec="(oddx4,sheet)  CHDx"
      integer:: Shape(2)
      integer,allocatable::  CHDIdx(:,:)

      integer:: FShape(2)
      integer,allocatable::  FCHDIdx(:,:)

      integer,pointer:: CHDnum(:,:)
      end       module modIndexing

      subroutine epSetIdx
      use modIndexing
      use modGetIndex
      implicit none
      integer i,j,  n, oddl
      integer Nlayer, NCHDx


      call epCountSubdTree(Spec, Shape, n, CHDnum)
      write(0,*) '1 num=', CHDnum(1,:)
      write(0,*) 'Shape=', Shape(:)

      oddl = index(Spec, "oddx4")
      if( oddl > 3 ) then
         write(0,*) 'Odd part contains',CHDnum(1,2), ' CHD,'
         write(0,*) 'they are located at the last part'
      else
         write(0,*) 'Odd part contains', CHDnum(1,1), ' CHD,'
         write(0,*) 'they are located at the fisrt part'
      endif
      Nlayer = Shape(1)
      NCHDx  = Shape(2)
      write(0,*) ' Nlayer=',Nlayer, ' NCHDx=',NCHDx
      allocate( CHDIdx(Nlayer, NCHDx) )

      call epGetIndex(Spec, Shape, CHDIdx, CHDnum)
      write(0,*) '2 num =', CHDnum(1,:)

      
      do i = 1, Nlayer
         do j = 1, NCHDx
            write(0,'(a, i4, a, i4, a, i6)')
     *      ' CHD # ',j, ' in layer ',i, ' has  comp#=',
     *      CHDIdx(i,j)
         enddo
      enddo

      forall(i=1:2) FShape(3-i) =Shape(i)
      allocate( FCHDIdx(FShape(1), FShape(2)))
      FCHDIdx =
     *      reshape(CHDIdx, FShape, order=(/2,1/))
      do j = 1, FShape(2)
         do i = 1, FShape(1)
            write(0,'(a, i4, a, i4, a, i6)')
     *      'F CHD # ',i, ' in layer ',j, ' has  comp#=',
     *           FCHDIdx(i,j)
         enddo
      enddo
      end
