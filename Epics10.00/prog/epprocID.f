      subroutine epprocIDins
!      when we find "id=..."  in a config line, this is colled
!      and apparent comp. # and id's are inserted with the 
!      struc name or subd name.   comp# is before expansion of 
!      subd.    
!
      use epidentify
      implicit none
#include  "Zep3Vec.h"
#include "Zcnfig.h"
!         now idis contain id list (comma separated items)
!     like   S1,layer1,trigger
      
      integer::idx
      character(len=snamelen)::name  
      integer::cpos,pos   ! comma pos in idis

      type(idblock),pointer:: block
      
      if(mode == 0) then
         idx = maxSubD + 1
      else
         idx = NsubD + 1
      endif
      name = Field(2)
      if( subDwithID(idx) == 0 ) then
         allocate( subdids(idx:idx) )
         block => subdids(idx)
         call epiniIdentify(block)
         subDwithID(idx) =1
      else
         if( associated( block) ) then
            nullify(block)
         endif
         block => subdids(idx)
      endif

      cpos = 1
      pos = 1
      do while (cpos  > 0) 
         cpos = index(idis(pos:), ",")
         if(cpos == 0 ) then
            call epinsertItem(
     *       block, trim(idis(pos:)), Det%nct, name )
         else
            call epinsertItem(
     *       block, idis(pos:cpos-1), Det%nct, name)
            pos = cpos + 1
         endif
      enddo
      end

      subroutine epprocIDincsubd
!
!      include and expand subdetector in the current sub 
!      detector or main detctor
!      
!
      use epidentify
      implicit none
#include  "Zep3Vec.h"
#include "Zcnfig.h"
!         now idis contain id list (comma separated items)
!     like   S1,layer1,trigger
      
      integer::idx
      character(len=snamelen)::name  

      character(len=idlen),pointer::ids(:)   ! id list 
      integer:: nids   !  # of id's   
      integer:: i, j
      type(idblock),pointer:: block
      
      if(mode == 0) then
         idx = maxSubD + 1
      else
         idx = NsubD + 1
      endif
      if( subDwithID(idx) > 0 ) then
!             subd  have id= or its children has
         if( associated(block) ) then
            nullify(block)
         endif
         block=>subdids(idx)
         if( associated(ids) ) then
            nullify(ids)
         endif
!         get the id's for the subd being expanded
         call epcn2idList(block, nnow, ids, nids)
!               change  the comp. no. to new one
         call epidChgNum(block, nnow, nnew)
         
         do i = 2,  SubD(subDx)%nct
            name = " "
            do j = 1, nids
               call epinsertItem(block, ids(j), nnew+i-1, name)
            enddo
         enddo
      endif
      end


