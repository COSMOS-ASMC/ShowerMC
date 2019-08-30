!      next one is   tentative sub for debug ; will  be detelet
      subroutine Lcompchk(id, no)
      implicit none
      character(*) id
      integer no
      if( no <= 0 ) then
         write(0, *)  " Lcomp (x) with x=", no, ' in ', id
!         stop
      endif
      end


!      collection of subs incluging #include
!    mixture of f90/ fixed format f77 with vax extension 
!  ****Be carefule*****
!     

      subroutine epLightCountDE
      use modepLight
      use modepLightPty
        !  assumed only called when Light > 0
        !  see if countde part has mn >=LinghtNoMin(=10), 
        !  if so, check if mn alredy appeared or not
        !      if first time, save mn and read propety file
        !
      implicit none
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"
#include "ZsepManager.h"
     
       integer::d, mn, B   !   CountDE: Bmnd is decomposed into element digit
       integer::i, j
       logical::exist
       logical::ok
       integer:: mdno
       integer:: mnCounter_tmp  ! from Akaike: 2010/12/7
       !                 Light Comp. count     Property count
       ! Light =11.            yes                yes
       !        12             yes                yes
       !        21             yes                yes (old no is NG)
       !        22             yes                yes                 
       !
       LightCompNo = 0 
       mnCounter = 0


       ok = .true.

       do i = 1, Det%nct
!                get d,m,B
          call epLightUnpackCountDE(Det%cmp(i)%CountDE, d, mn, B)

          if( mn < LightNoMin )  then
             Det%cmp(i)%LightCompNo = 0 ! this comp. is not related to light
          else
             if( LightCompNo >= maxLightComp )  then
                ! too many light components
                write(0, *)
     1             ' too many components for light related comp>',      
     2             maxLightComp
                write(0, *) 
     &            ' may be you must increase maxLightComp in ',
     &            ' ZlightMaxDef%h'
                stop 0000
             endif
             LightCompNo = LightCompNo + 1

             Lcomp(LightCompNo)%mn = mn
             Det%cmp(i)%LightCompNo = LightCompNo
             Lcomp(LightCompNo)%compno = i

             exist = .false.    ! for mn check
             do j = 1, mnCounter
                if( mn == mnContainer(j) ) then
                   exist=.true.   ! same as existing one
                   exit
                endif
             enddo
             if(exist) then
!                 mn already in mnContainer
!cc                if(Light /= 21) then
                   ! need ray tracing
                   Lcomp(LightCompNo)%comInfoNo = j  ! com info is here
                                                   ! i.e, comInfo(j)
!                   Lcomp(LightCompNo)%mn = mn
!cc                else
!cc                   write(0,*)
!cc     &               ' strange Light=21==>property file not needed'
!cc                   write(0,*) 
!cc     &              ' so we should not come here; in  epLightCountDE'
!cc                   ok=.false.
!cc                endif
             else
          !          new mn
!cc                if(Light /= 21 ) then
                   if(mnCounter >= maxProperties) then
                      write(0,*)
     &                 ' too many diff. mn in CountDE>',
     &                   maxProperties
                      write(0,*) ' see maxProperties in epLightMaDef%h'
                      stop 5432
                   endif
                   mnCounter = mnCounter + 1
                   mnContainer(mnCounter) = mn

                   Lcomp(LightCompNo)%comInfoNo = mnCounter
!                   Lcomp(LightCompNo)%mn = mn
  !               if Light /=21, we  must read component and media 
            !    property file  for light generation and tracing
                   call epLightParamRead0(LightDir, mn,  
     &              comInfo( mnCounter ) )
                   if( comInfo(mnCounter)%refracIndex(1) == 0.) then
               !     this is spcial. if refracindex is not given, use
               !     madia value; 
                      mdno = Det%Cn2media(i)
                      comInfo(mnCounter)%refracIndex = Media(mdno)%n
                   endif
            !      read wave length distribution file etc, if spcified

                   call epLightReadAFFile(LightDir, Lcomp(LightCompNo),
     &                   comInfo(mnCounter) )
!cc                endif
             endif
               !  compute some values and store Lcomp
               !  we should re-fix the two index values here 
               !  for known mn case

             LightCompNo =  Det%cmp(i)%LightCompNo 

!             mnCounter =  Lcomp(LightCompNo)%comInfoNo 
             mnCounter_tmp =  Lcomp(LightCompNo)%comInfoNo   ! from Akaike
             call epLightManipInp(i,
!     &           comInfo(mnCounter), Lcomp(LightCompNo) )
     &           comInfo(mnCounter_tmp), Lcomp(LightCompNo) )  ! //
          endif
       enddo


       if(LightCompNo == 0 ) then
          write(0,*) 
     &  ' Light > 0 but no component with mn of CountDE> 0'
          stop 9999
       else
          write(0,*) '# of Light components is ', LightCompNo
          write(0,*) '# of unique property files is ', mnCounter
       endif
       if(.not. ok) then
          write(0,*) ' There is some wrong spcification in '
          write(0,*) ' CountDE part of config '
          stop 9999
       endif
       end  

!           component # to Lcomp #
      subroutine epLightCn2LcompNo(cmpNo, LcmpNo)
      implicit none

#include "ZepTrackv.h"
#include "Zcnfig.h"
      integer cmpNo  ! input comp. #
      integer Lcmpno ! output. cn-->LightCompNO 
                    !  If not Light component, 0
      LcmpNo = Det%cmp(cmpNo)%LightCompNo
      end
      subroutine epLightLcompNo2Cn(LcmpNo, cmpNo)
      implicit none
#include "ZepTrackv.h"
#include "Zcnfig.h"
      integer,intent(in):: Lcmpno ! LightCompNO 
      integer,intent(out)::cmpNo  !  comp. #; if Lcmpno is fake, cmpNo=0
      logical,save::first=.true.
      integer,allocatable,save::Lc2c(:)
      integer ic, ilc
      if( first ) then
         allocate( Lc2c(Det%nct) )
         first = .false.
         Lc2c(:) = 0 
         do ic = 1, Det%nct
            ilc = Det%cmp(ic)%LightCompNo
            if(ilc > 0 ) then
               Lc2c(ilc) = ic
            endif
         enddo
      endif
      if(LcmpNo > 0 ) then
         cmpNo = Lc2c(LcmpNo)
      else
         cmpNo = 0
      endif
      end

      
!        set  countDE within program after reading config data.
!        may be called by the user from uiav.
!        If this is called, epLightCountDE( user ) must be
!        called ; as to user see epLightCounDE  
      subroutine epLightResetCountDE(cmpNo, cde)
      implicit none

#include "ZepTrackv.h"
#include "Zcnfig.h"
      integer,intent(in):: cmpNo  !  comp. #
      integer,intent(in):: cde    ! . countDE to be set
      if(cmpNo >= 1 .and. cmpNo <= Det%nct) then
         Det%cmp(cmpNo)%CountDE = cde
      else
         write(0,*) 'Warning from epLightResetCountDE:'
         write(0,*) ' specified comp.# ', cmpNo, ' non exsistent'
      endif
      
      end

      subroutine epLightMn2cominfoIdx(mn, idex)
      use modepLight
      implicit none
!     "mn" in  countDE  part of a scintillator specifies
!      properies of the scintillator is in file Lightmn.dat
!      in "LightDir".  (If mn > 50, not scintillator but
!      is a sensor.)  The properties may be common to 
!      a number of scintillators and such common info. is
!      stored in  comInfo( idex )%...  where  
!      mnConainer(idex)= mn

      integer,intent(in):: mn 
      integer,intent(out):: idex   ! will be 0 if mn is non existent

      integer i

      do i = 1, mnCounter 
         if( mnContainer(i) == mn ) then
            idex = i
            exit
         endif
         idex = 0
      enddo
      end


