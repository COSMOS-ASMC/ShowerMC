      module modcountStruc
      implicit none
#include  "ZepMaxdef.h"
!               max number of diff. structur names (such as
!           box octagon octagon_y  et. 
!           5 is for not 'new' one (box, pipe, prism, cyl, sphere).
!           3 is for _y _z 
      integer,parameter::nname=(MAX_NEW_STRUC+5)*3
      character(MAX_STRUCCHR)::namelist(nname)
      integer,save:: counter(nname) ! counter for each struc
      integer,save:: indx(nname) ! used for sorting
      integer,save:: used    ! counter of diff. struc
      logical,save:: error   ! error flag for too many diff. struc
      end   module modcountStruc

      subroutine epcountStruc
      use modcountStruc
      implicit none
#include  "ZepMaxdef.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"

      integer  i, j
      counter(:) = 0
      indx(:) = 0

      namelist(1) = 'box'   ! box is so many ; to reduce 
               ! computation time, put it at the first place
      used = 1   !  but counter  is still 0
!      Det.nworld = 0    ! has been set by recount
      error = .false.
      do  i = 1, Det%nct  ! for all components
         Det%cmp(i)%cn = i  ! renumber
         do j = 1, used     ! check already appeared ?
            if( namelist(j) ==  Det%cmp(i)%struc ) then
               counter(j) = counter(j) + 1  ! count up
               goto 100  ! see next comp.
            endif
         enddo
         ! new name appeared
         if(used >= nname) then
            write(0,*) ' too many diff.  structure names '
            write(0,*) 
     *   ' you have to increase MAX_STRUCCHR in ZepMaxdef%h'
            write(0,*) ' or nname in eprecount sub in eprcnf%f'
            error = .true.
            exit
         endif
         used = used + 1
         namelist(used) = Det%cmp(i)%struc
         counter(used) = 1
 100     continue
      enddo
        !  sort by name
      if(counter(1) == 0 ) then
         ! no box appeared; move the list 
         do i = 2, used
            namelist(i-1) = namelist(i) 
            counter(i-1) = counter(i)
         enddo
         used = used - 1
      endif
!        asecnding sort by name  
      call kqsortc(namelist(1), indx, used)
      end     subroutine epcountStruc
      
      subroutine epwriteStruc
!        show list of each struc(vol-shape) and # of use
      use   modcountStruc
      implicit none
#include  "ZepMaxdef.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer :: i

      write(0,*) ' vol(shape)-name : # of appearance'
      do i = 1, used
         write(0, '(4x,a12,3x,i7)') 
     *         namelist(indx(i)), counter( indx(i) )
      enddo
      if(error) stop
!          this should have been done in eprecount
!      if(index(trim(Det.cmp(Det.nct).struc), '_w') .ne. 0) then
!         Det.nworld = 1
!      endif
      
      end   subroutine epwriteStruc

      subroutine epqNoOfThisStruc(name, num)
!       ask how many times this 'name' is used
      use   modcountStruc
      implicit none
      character(len=*),intent(in):: name ! search # of stuctures with                               ! this name
      integer,intent(out):: num  ! >=0.  # of appearance of such 
                                 ! volshapes
                                 
      integer::i
      num = 0
      do i = 1, used
         if(namelist(i) == name ) then
            num = counter(i)
            return
         endif
      enddo
      end  subroutine epqNoOfThisStruc
