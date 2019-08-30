#include "ZepicsBD.h"
#include "ZcosmosBD.h"
      module modtree2idx
      integer::n
      integer:: shape(7)
      integer,pointer:: num(:,:)
      integer,allocatable:: Idx1(:)
      integer,allocatable:: Idx2(:,:)
      integer,allocatable:: Idx3(:,:,:)
      integer,allocatable:: Idx4(:,:,:,:)
      integer::     j1, j2, j3, j4
      character(len=16):: target
      character(len=100):: spec
      character(len=5):: nth(4)=(/"1st","2nd","3rd","4th"/) 
      end       module modtree2idx
      program tree2index
      use modGetIndex, nf=>nitem
      use modtree2idx
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
      character*100  dsn1
      integer:: i
      MediaDir(1) = '$EPICSTOP/Data/Media'
      write(0,*) 'Enter  config file path'
      read(*,'(a)') dsn1
      write(0,*) 'Enter, subd hierarchy: eg:  sheetX  SciFi'
      read(*,'(a)') spec
      write(0,*) ' reading config'
      call eprcnf(dsn1)

      shape(:) = 0
      call epCountSubdTree(spec, shape, n, num)
      if(n == 1 ) then
         allocate( Idx1(shape(1) ) )
         call epGetIndex(spec, shape, Idx1, num)
      elseif(n == 2 ) then
         allocate( Idx2(shape(1), shape(2) ))
         call epGetIndex(spec, shape, Idx2, num)
      elseif( n == 3) then
         allocate( Idx3(shape(1), shape(2), shape(3) ))
         call epGetIndex(spec, shape, Idx3, num)
      elseif( n == 4) then
         allocate( Idx4(shape(1), shape(2), shape(3), shape(4) ))
         call epGetIndex(spec, shape, Idx4, num)
      else
         write(0,*) 'too many hierarchy;n=',n, ' must be <= 4'
         stop
      endif
      write(*,*) 'spec=', trim(spec)
      write(*,*) 'shape info ', shape(1:n)
      write(*,*) 'num info' 
      write(*,*) 'i  num(i, 1)   num (i, 2)'        
      do i = 1, ntbr
         write(*,'(i2, 2i9)') i, num(i,1:OddSubd)
      enddo

      target = splitSpec(nf)
      if( ntbr == 1 ) then
         if( OddSubd == 1 ) then
            call tree2idx11h
         else
            call tree2idx12h
         endif
      else
         if( OddSubd == 1 ) then
            call tree2idx21h
         else
            call tree2idx22h
         endif
      endif
!           print index table
      select case(n)
      case(1)
         call tree2idxCase1
      case(2)
         call tree2idxCase2
      case(3)
         call tree2idxCase3
      case(4)
         call tree2idxCase4
      case default
         write(0,*) ' strange error'
         stop
      end select
      end  program tree2index

      subroutine tree2idx11h
!          header for ntbr=1 OddSubd=1 case
      use modGetIndex
      use modtree2idx
      implicit none
        ! plain case ; spec is "A B C" type
      write(*,'(a,a,a)') 
     *     '1st index is for ',trim( topbranch(1) )
      write(*,'(a,a,a,i5,a)')
     *      "# of ",trim( topbranch(1)), " is ", shape(1),
     *      " [=shape(1)]"
      write(*,'(a,a,a,i4,a,a,a)')
     *      "In each ",trim( topbranch(1)),", ", shape(n),
     *      " [=shape(n)] ",trim(target), " are contained "
      write(*,'(a,a,a,i6,a)')          
     *      "total # of ",trim(target)," is ",
     *     product( shape(1:n) )," [=product( shape(1:n) )]"
      end subroutine tree2idx11h

      subroutine tree2idx12h
!          header for ntbr=1 OddSubd=2 case
        !  spec is "A  (B1,B2) C" type nf>2
        ! or
        !          "(B1,B2) C" type nf==2
      use modGetIndex, nf=>nitem
      use modtree2idx
      implicit none

      if(nf==2) then
         call tree2idx12h_nf2
      else
         call tree2idx12h_nf3
      endif
      end       subroutine tree2idx12h
      subroutine tree2idx12h_nf2
!          header for ntbr=1 OddSubd=2 case
        !          "(B1,B2) C" type nf==2
      use modGetIndex, nf=>nitem
      use modtree2idx
      implicit none
      integer:: i

      write(*,'(a,i2,a)' )
     *  '1st idx varies from 1 to ',shape(1),' [=shape(1)]'
      write(*, '(a)') 'It may be layer # or something ...'
      write(*,'(a,a,a,i4,a)') 
     *  '# of ',trim(target),' in each "layer" is ',
     *   shape(2),' [=shape(2)]'
      write(*,'(a)') 'ingredient is '
      write(*, '(a,a)') 'subdname     # of ',trim(target)
      write(*,'(a,i4,a)')   trim(OddSubdName(1)), 
     *          num(1,1), " [=num(1,1)]"
      write(*,'(a,i4,a)')   trim(OddSubdName(2)), 
     *          num(1,2), " [=num(1,2)]"

      write(*,'(a,a,a,i6,a)')          
     *      "total # of ",trim(target)," is ",
     *       (num(1,1)+num(1,2))*shape(1),
     *     " [=(num(1,1)+num(1,2))*shape(1)]"

      end subroutine tree2idx12h_nf2
      subroutine tree2idx12h_nf3
!          header for ntbr=1 OddSubd=2 case
        !          "(B1,B2) C" type nf==2
      use modGetIndex, nf=>nitem
      use modtree2idx
      implicit none
      integer:: i
      end subroutine tree2idx12h_nf3

      subroutine tree2idx21hxx ! not used. see  21h below
!          header for ntbr>1 OddSubd=1 case
        !  spec is "|A1,A2..|  B C" type
      use modGetIndex, nf=>nitem
      use modtree2idx
      implicit none
      integer:: i
      write(*, '(a)') '1st index    subdname        #  '
      do i = 1, ntbr 
         write(*,'(i5, a,a,i4,a)')  i, " ", trim(topbranch(i)), 
     *           shape(2), " [=shape(2)]"
      enddo
      write(*,'(a,i4,a,a,a)')
     *      "In each branch ", shape(n),
     *      " [=shape(n)] ",trim(target), " are contained "
      write(*,'(a,a,a,i6,a)')          
     *      "total # of ",trim(target)," is ",
     *     product( shape(1:n) )," [=product( shape(1:n) )]"

      end subroutine tree2idx21hxx

      subroutine tree2idx21h
!          header for ntbr>1  OddSubd=1 case
        !  spec is "|A1,A2| B C" type nf>2
        !          "!B1,B2| C" type nf==2
      use modGetIndex, nf=>nitem
      use modtree2idx
      implicit none
      integer:: i
      write(*,'(a,a)')
     *      " 1st idx   subdname   # of ",trim(target)
      if( shape(1) /= ntbr) then
         write(0,*) ' ntbr=',ntbr, ' /= shape(1)=',shape(1)
         stop
      endif

      do i = 1, ntbr
         write(*,'(i4, 7x, a, i8, a,i1,a)' ) 
     *      i, trim(topbranch(i)), num(i,1), " [=num(",i,",1)]"
      enddo
      end    subroutine tree2idx21h

      subroutine tree2idx22h
!          header for ntbr=2 OddSubd=1 case
        !  spec is "|A1,A2| B C" type nf>2
        !          "!B1,B2| C" type nf==2
      use modGetIndex, nf=>nitem
      use modtree2idx
      implicit none
      integer:: i
      write(*,'(a,a)')
     *      " 1st idx   subdname   # of ",trim(target)
      do i = 1, ntbr
         write(*,'(i4, 2x,a, i4,a)' ) 
     *      i, trim(topbranch(i)), num(i,1), " [=num(i,1)]"
      enddo
      end    subroutine tree2idx22h


      subroutine tree2idxCase1
        ! index and  Cn list for n=1
      use modGetIndex
      use modtree2idx
      implicit none
      integer:: i
      
      write(*,'(a,a)') "    index     ", trim(target)
      write(*,'(3x,a,a)' )
     *      (nth(i), i=1,n), " Cn  # "
      do j1 = 1, shape(1)
         if( OddSubd >1) then
            write(*,*) "=====================",trim(OddSubdName(1))
         endif
         write(*,'( i5, i8)')
     *           j1,  Idx1(j1)
      enddo
      end subroutine tree2idxCase1

      subroutine tree2idxCase2
        ! index and  Cn list for n=2
      use modGetIndex
      use modtree2idx
      implicit none
      integer:: i
      
      write(*,'(a,a)') "    index     ", trim(target)
      write(*,'(3x,a,a,a)' )
     *      (nth(i), i=1,n), " Cn  # "
      do j1 = 1, shape(1)
         if( OddSubd >1) then
            write(*,*) "=====================",trim(OddSubdName(1))
         endif
         do j2 = 1, shape(2)
            if( j2 == num(1,1)+1 .and. OddSubd>1) then
               write(*,*) '------------------',trim(OddSubdName(2))
            endif
            write(*,'( i5, i5, i8)')
     *           j1, j2,  Idx2(j1,j2)
         enddo

      enddo
      end subroutine tree2idxCase2

      subroutine tree2idxCase3
        ! index and  Cn list for n=3
      use modGetIndex
      use modtree2idx
      implicit none
      if( ntbr == 1) then 
         call tree2idxCase3tb1
      else
         call tree2idxCase3tb2
      endif
      end  subroutine tree2idxCase3

      subroutine tree2idxCase3tb1
        ! index and  Cn list for n=3 ntbr=1
      use modGetIndex
      use modtree2idx
      implicit none
      integer:: i

      write(*,'(a,a,a)') "    index     ", trim(target)
      write(*,'(3x,a,a,a,a)' )
     *      (nth(i), i=1,n), " Cn  # "
      
      do j1 = 1, shape(1)
         write(*,'(a)') '=================',trim(splitSpec(1))
         do j2 = 1, shape(2)
            if( OddSubd >1) then
               write(*,*)
     *           "-------------------",trim(OddSubdName(1))
            else
               write(*,*)
     *           "-------------------",trim(splitSpec(2)) 
            endif
            do j3= 1, shape(3)
               if( j3 == num(1,1)+1 .and. OddSubd>1) then
                  write(*,*)
     *             '........',trim(OddSubdName(2))
               else
                  write(*,*)
     *             '........',trim(splitSpec(3))
               endif
               write(*,'( i5, i5, i5, i8)')
     *              j1, j2, j3,  Idx3(j1,j2,j3)
            enddo
         enddo
      enddo 
      end subroutine tree2idxCase3tb1

      subroutine tree2idxCase3tb2
        ! index and  Cn list for n=3 ntbr>1
      use modGetIndex
      use modtree2idx
      implicit none
      integer::i

      write(*,'(a,a)') "      index       ", trim(target)
      write(*,'(3x,a,a,a,a)' )
     *      (nth(i), i=1,n), " Cn  # "
      
      do j1 = 1, shape(1)
         write(*,'(a,a)') '===',trim(topbranch(j1))
         do j2 = 1, shape(2)
            write(*,'(a)') '-------"layer"' 
!            do j3= 1, shape(3)
            do j3= 1, num(j1,1)
               if( j3 == 1 ) then
                  write(*,'(a,a)')
     *           "............",trim(splitSpec(2)) 
               endif                  
               write(*,'( i5, i5, i5, i8)')
     *              j1, j2, j3,  Idx3(j1,j2,j3)
            enddo
         enddo
      enddo 
      end subroutine tree2idxCase3tb2



      subroutine tree2idxCase4
      use modGetIndex
      use modtree2idx
      implicit none
      integer:: i 
      write(*,'(a,a)') "      index          ", trim(target)
      write(*,'(3x,a,a,a,a,a)' )
     *      (nth(i), i=1,n), " Cn  # "

      do j1 = 1, shape(1)
         do j2 = 1, shape(2)
            do j3= 1, shape(3)
               do j4 = 1, shape(4)
                  write(*,'( i5, i5, i5, i5, i8)')
     *              j1, j2, j3, j4,  Idx4(j1,j2,j3,j4)
               enddo
            enddo
         enddo
      enddo
      end  subroutine tree2idxCase4
