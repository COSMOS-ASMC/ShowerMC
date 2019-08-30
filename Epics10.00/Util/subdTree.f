      module modsubdTree
      implicit none
      integer,save::namelen = 0      
      integer,save::option= 0
      end module modsubdTree
      subroutine epSubdTree(out)
      use modsubdTree
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!         
       integer,intent(in):: out  !   output logical device no.
       integer::i
       write(out,*) 'L---- subdetector list ----'
       do i = 1, NsubD
          write(out,*) i, subdName(i)
          namelen = max(namelen, len(trim(subdName(i))))
       enddo
       write(out,*) ' '
       write(out,*) 'H---- subdetector hierarchy ----'
       call subdTree0(out, 0, Det%nct)

       end

       recursive subroutine subdTree0(out,mother, i)
       use modsubdTree
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
       integer,intent(in)::out ! logical output dev #
       integer,intent(in)::i  ! comp. #
       integer,intent(in)::mother ! subd index of the 
                            !mother component of comp. # i
                ! if  i is  world, 0

       integer:: depth
       integer:: j, m, k, ii, l
       integer:: indexs       
       character(17):: sname
       character(100)::line
       integer:: msubdi, flag, ll

       depth=Det%cmp(i)%level
       j = Det%cmp(i)%NMatreska
       if( j .gt. 0 .or. option==1 )  then
          call epqSubdName2(mother, i, sname, flag)
          line = ' '
          ll =len( trim(sname) ) + 1
          if( flag == 0 ) then
             sname=trim(sname)//"*"
          else
             sname=trim(sname)//" "
          endif
          write(line(3*depth+1:),
     *       '("|",i2,1x, a,i7,1x, a,a)') 
     *       depth,  sname(1:max(namelen,ll) ),
     *          i,  Det%cmp(i)%struc,
     *         Det%cmp(i)%matter
          write(out,*) trim(line)
          call epqSubdIdx(i, msubdi)
          do m = j, 1, -1
                 ! for partially one, CnArea.... < 0
             ii =abs( CnArea( Det%cmp(i)%ContainsR+m ) )
             call subdTree0(out, msubdi,  ii)
          enddo
       endif


       end subroutine subdTree0


