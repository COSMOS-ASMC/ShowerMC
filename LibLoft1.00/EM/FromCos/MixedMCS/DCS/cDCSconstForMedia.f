       subroutine cDCSConstForMedia(md, pm, aDCS)
!  use modXSecMedia
       use modTPXS
       use modDCS
       implicit none
#include "Zmedia.h"
       
       type(epmedia),intent(in):: md
       integer,intent(in):: pm   ! >0 for e+, < 0 for e-
       
       type(DCSconst),intent(inout):: aDCS  ! DCSconst for media.
       
       type(DCSconst),pointer::zDCS(:) ! for each Z in media
                                 ! allocated here and deallocated 

       integer:: n  ! # of elements in the medium
       integer:: nmuc, nec  

       integer,parameter::io=11  ! temporaray disk logcial dev. #
                           ! should be the same as TempDeV
       integer:: i, j, k,  Z, icon, nE
       real(8):: temp

!       to fix KEele, nEele, nEpos in modTPXS; now not needed
!    if( ( pm < 0 .and. nEneg == 0 ) .or. &
!        ( pm > 0 .and. nEpos == 0 )  ) then
!       to fix KEele

         call cPreReadDCS(io)
!    endif

         n = md%noOfElem
         allocate(zDCS(n))

         nE = nEneg
         if( pm > 0 ) nE = nEpos
         do i = 1, n
            Z = md%elem(i)%Z
                  ! create e-/e+ DCS file
            do nec = 1, nE
               call cDCS_Z2file(Z, pm, nec, filename)
               call copenf(io, filename, icon) 
               if(icon /= 0 ) then
                  write(0,*) ' file=',trim(filename),
     *             ' cannot be opened '
                  stop
               endif
               call cReadDCS(io,pm, nec, zDCS(i))  ! for each Z and E, set dcs
            enddo
         enddo
               

         do nec = 1, nE
            do nmuc = 1, nmu
               temp = 0.
               do k = 1, n
                  temp = temp +md%No(k) * zDCS(k)%dcs(nmuc, nec)
               enddo
               aDCS%dcs(nmuc, nec) =temp
            enddo
         enddo
         deallocate(zDCS)
         write(0,*) ' entering cPrepIntpDCs'
         call cPrepIntpDCS(aDCS)
      
       end subroutine cDCSConstForMedia
