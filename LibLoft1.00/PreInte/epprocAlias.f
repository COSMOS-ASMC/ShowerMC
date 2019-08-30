#include  "ZepMaxdef.h"
      module modAlias
      implicit none
      integer,parameter::maxleng = MAX_MEDIANAMELENG
      integer,save:: NoOfAlias = 0
      integer,parameter:: MaxNoOfAlias = 20
      character(maxleng),save::
     *  aliasName(MaxNoOfAlias) = " "
      character(maxleng),save::
     *  trueName(MaxNoOfAlias) = " "
      end module modAlias

      subroutine epProcAlias(Field)
      use modAlias
      implicit none
      character(*),intent(in)::Field(3) ! fieldAsItis
      integer:: i

!        check if alias is already used.

      do i = 1,  NoOfAlias 
         if( trim( Field(2) ) == trim( aliasName(i) )) then
            if( trim( truename(i) ) == trim( Field(3)) ) then
               write(0,*) 
     *         "******************************",
     *         trim( Field(1))," ", trim( Field(2))," ",
     *         trim( Field(3)), " already defined"
               write(0,*) " so you may omit it "
               return !!!!!!!!!!!!!!
            else
               write(0,*)
     *        'alias name=',trim( Field(2)), ' already used for ',
     *        trim( truename(i)), ' so cannot be used for ',
     *        trim(  Field(3))
               stop
            endif
         endif
      enddo
         
      if( NoOfAlias < MaxNoOfAlias ) then
         NoOfAlias =  NoOfAlias + 1
      else
         write(0,*) 'too many  aliases for the media name'
         write(0,*)
     *     ' Enlarge MaxNoOfAlias in modAlias in eprcnf%f '
         write(0,*) ' it is now ',MaxNoOfAlias
         stop
      endif
!         check of trueName will be done after all
!         config data is read.
      aliasName(NoOfAlias) = Field(2)
      trueName(NoOfAlias) = Field(3)
      end       subroutine epProcAlias
      
      subroutine epSeeIfAlias(name)  
!       see if name is alias. if so, it is replace by
!       the true name
      use modAlias
      implicit none
      character(*),intent(inout):: name
      integer::i
      do i = 1, NoOfAlias
         if( trim( aliasName(i)) == trim( name) ) then
            name = trim( trueName(i) )
            exit
         endif
      enddo
      end   subroutine epSeeIfAlias

      subroutine epwriteAliasMedia
      use modAlias
      implicit none
!          print media info
! #include "ZepTrackv.h"
! #include "Zcnfig.h"
! #include "ZepManager.h"

      integer:: i
      if( NoOfAlias > 0 ) then
         write(0,*) 'alias  of  media' 
      endif
      do i = 1, NoOfAlias
         write(0,*) aliasName(i),"   ", trueName(i)
      enddo
      end     subroutine epwriteAliasMedia
      
