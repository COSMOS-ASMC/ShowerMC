!       For testing odd detector attachment 
!          configB3 must be used;  to be specified in FirstInput
      module modIndexing
      character(len=*),parameter:: 
     *  Spec="|sheet,odddummy| CHDx" 
      integer:: Shape(3)
      integer,allocatable::  CHDIdx(:,:,:)
      integer,pointer:: CHDnum(:,:)

      integer:: FShape(3)
      integer,allocatable::  FCHDIdx(:,:,:)
      end       module modIndexing

      subroutine epSetIdx
      use modIndexing
      use modGetIndex
      implicit none
      integer i,j,  n
      integer Nlayer, NCHDx

      call epCountSubdTree(Spec, Shape, n, CHDnum)
      write(0,*) ' num(1,1), num(2,1)=', CHDnum(1,:), CHDnum(2,:)
      write(0,*) 'Shape=', Shape(:)

      Nlayer = Shape(2)
      NCHDx  = Shape(3)
      write(0,*) ' Nlayer=',Nlayer, ' NCHDx=',NCHDx
      allocate( CHDIdx(Shape(1),Nlayer, NCHDx) )

      call epGetIndex(Spec, Shape, CHDIdx, CHDnum)
      
      do i = 1, Nlayer
         do j = 1, NCHDx
            write(0,'(a, i4, a, i4, a, i6)')
     *      ' CHD # sheet ',j, ' in layer ',i, ' has  comp#=',
     *      CHDIdx(1, i,j)
         enddo
      enddo

      do i = 1, Nlayer
         do j = 1, NCHDx
            write(0,'(a, i4, a, i4, a, i6)')
     *      ' CHD # odddummy ',j, ' in layer ',i, ' has  comp#=',
     *      CHDIdx(2, i,j)
         enddo
      enddo


      forall(i=1:3) FShape(4-i) =Shape(i)
      allocate( FCHDIdx( FShape(1), FShape(2), Fshape(3)) )
      FCHDIdx =
     *      reshape(CHDIdx, FShape, order=(/3,2,1/))
      do j = 1, FShape(2)
         do i = 1, FShape(1)
            write(0,'(a, i4, a, i4, a, i6)')
     *      'F CHD # sheet ',i, ' in layer ',j, ' has  comp#=',
     *           FCHDIdx(i,j,1)
         enddo
      enddo
      end
