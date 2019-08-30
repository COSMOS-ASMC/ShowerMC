      subroutine epqevn(nev)
      implicit none
#include  "ZepManager.h"
!            inquire  current event number created in this run
       integer nev
       nev = Nevrun
       end
      subroutine epqCn2Media(compn, mediax)
  !  component number to media
      implicit none
#include "ZepTrackv.h"
#include "Zcnfig.h"
      integer,intent(in):: compn !  component number
       type(epmedia)::  mediax  ! output.

      call epchkcmpn('epqCn2MediaIdx', compn)
      mediax = Media( Det%Cn2media(compn) )
      end
      subroutine epqSubdIdx(compn, idx)
  !  component number to subdetector index
      implicit none
#include "Zep3Vec.h" 
#include "Zcnfig.h"
      integer,intent(in):: compn !  component number
      integer,intent(out):: idx  ! subdetector index to which
                     ! the component belongs (if compn 
                     ! is a simple comp.)
                     ! If it is a subdetector, idx becomes 
                     ! its subd index. 
                     ! if compn is invalid, -1
                     ! if compn does not belong to a subd, 0 
      if(compn >=1 .and. compn <= Det%nct) then
         idx=Det%cmp(compn)%subdidx
      else
         idx = -1
      endif 
      end   subroutine epqSubdIdx
      subroutine epqSubdIdx2(compn, idx)
  !  component number to subdetector index
      implicit none
#include "Zep3Vec.h" 
#include "Zcnfig.h"
      integer,intent(in):: compn !  component number
      integer,intent(out):: idx  !  if compn is a simple
               !          object, 0
               ! If it is a subdetector, or simple object
               ! containg children, idx becomes 
               ! its subd index. 
               ! if compn is invalid, -1

      if(compn >=1 .and. compn <= Det%nct) then
         if( Det%cmp(compn)%NMatreska > 0 ) then
            idx=Det%cmp(compn)%subdidx
         else
            idx = 0
         endif
      else
         idx = -1
      endif 
      end   subroutine epqSubdIdx2

      subroutine epqSubdName(compn, subdetectorName)
  !    get subdetector name of given component number
  !  If the component is a simple object, 
  !  it becomes  the subdetecor name which contain the object
  !  If it is not contained, ' '

  !  If the component is a subdetector, it becomes that 
  !  subdetector name. (cf epqMotherName)

      implicit none
#include "Zep3Vec.h" 
#include "Zcnfig.h"
      integer,intent(in):: compn ! component #
      character(*),intent(out)::  subdetectorName ! must be >=16 char.
                         ! subdname; may be blank, if compn 
                         ! does not belong to ' '  

      integer:: idx
      call epqSubdIdx(compn, idx)
      if(idx > 0 ) then
         subdetectorName = SubDname(idx)
      elseif( index(Det%cmp(compn)%struc, "_w") > 0 )then
         subdetectorName = 'world '
      else
         subdetectorName = ' '
      endif
      end  subroutine epqSubdName

      subroutine epqSubdName2(mother, compn, Name, flag)
  !   get subdetector name of given component number
  !  If the component is a simple object, 
  !  it becomes  the structure (box etc) name (cf epqSubdName)

  !  If the component is a subdetector, it becomes that 
  !  subdetector name. (cf epqMotherName)

      implicit none
#include "Zep3Vec.h" 
#include "Zcnfig.h"
      integer,intent(in):: mother ! subd index of
              !  mother comp# of compn
              !  if 0, compn is in the world.
      integer,intent(in):: compn ! component #
      character(*),intent(out)::  Name ! must be >=16 char
      integer,intent(out):: flag  ! if compn is simple 
             ! object 0, if subdetector, 1
      integer:: idx
      call epqSubdIdx2(compn, idx)
      if(idx > 0 ) then
         if( idx == mother) then
            ! should be a simple object
            Name  =Det%cmp(compn)%struc
            flag = 0
         else
            flag = 1
            Name = SubDname(idx)
         endif
      elseif( index(Det%cmp(compn)%struc, "_w") > 0 )then
         Name = 'world '
         flag = 1
      else
         Name = Det%cmp(compn)%struc
         flag = 0
      endif

      end  subroutine epqSubdName2

      subroutine epqMotherName(compn, mindex, motherName)
  !    get subdetector name which contains the given component number
  !    If it is not contained, ' '
! *** If you need to use this one during simulation frequently,
!  it's better to make hash table at the init time of all the events 
      implicit none
#include "Zep3Vec.h" 
#include "Zcnfig.h"
      integer,intent(in):: compn ! component #
      integer,intent(out):: mindex ! mother's index for subdName
      character(*),intent(out)::  motherName ! must be >=16 char.
                         ! subdname; may be blank, if compn 
                         ! does not belong to ' '  

      integer:: idx, i, j, m, ii
      do i = Det%nct, 1, -1
         j = Det%cmp(i)%NMatreska
         if( j >  0 ) then
            do m = j, 1, -1
               ii =abs( CnArea( Det%cmp(i)%ContainsR+m ) )
               if( ii == compn )  then
                  call epqSubdIdx(i, mindex)
                  return  !*************
               endif
            enddo
         endif
      enddo
      
      end  subroutine epqMotherName
      
      subroutine epqCompsInSubD(subdidx, ncomps)
  !    inquire the number of components in a given
  !    subdetector specified by the index
      implicit none
#include "Zep3Vec.h" 
#include "Zcnfig.h"
      integer,intent(in):: subdidx ! subdetecor index
                                ! if <=0, main is assumed

      integer,intent(out):: ncomps ! # of comps. in the subd
                                ! if subdidx  > max SubD, 0 
      if(subdidx == 0 ) then
         ncomps = Det%nct
      elseif( subdidx > NsubD) then
         ncomps = 0
      elseif( subdidx < 0 ) then
         ncomps = 0
      else
         ncomps = SubD(subdidx)%nct
      endif
      end    subroutine epqCompsInSubD

      subroutine epqNoOfMedia(nmedia)
  !    total # of media used in the current job
      implicit none
#include "ZepTrackv.h"
      integer,intent(out):: nmedia ! # of media
      nmedia = NoOfMedia
      end
      subroutine epqTotal1ryE(sumKE)
  !    total K.E of the 1ry pariticles of the current event
      implicit none
#include "ZepTrackv.h"
      real(8),intent(out):: sumKE
      sumKE = Total1ryE
      end
      subroutine epqCn2MediaName(compn, medianame)
      implicit none
#include "ZepTrackv.h"
#include "Zcnfig.h"
      integer,intent(in):: compn 
      character(*),intent(out):: medianame  ! >= 8 cha

      integer:: mediaidx ! if compn wrong, 0 else 
                    ! its media index is obtained
      call epqCn2MediaIdx(compn, mediaidx)
      if( mediaidx > 0 ) then
         medianame = Media(mediaidx)%name
      else
         medianame = 'xxxxxxx'
      endif
      end
      subroutine epqCn2AliasMediaName(compn, aliasname, name)
      implicit none
#include "Zep3Vec.h"      
#include "Zcnfig.h"      
      integer,intent(in):: compn
      character(*),intent(out):: aliasname ! >= 8 cha.  if 
          ! alias is not used, same as name.
      character(*),intent(out):: name ! >= 8 cha . true media name
      call epqCn2MediaName(compn, name)
      aliasname = Det%cmp(compn)%matter
      end  subroutine epqCn2AliasMediaName

      subroutine epqCn2MediaIdx(compn, mediaidx)
!         component # to media index
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
      integer,intent(in):: compn 
      integer,intent(out):: mediaidx ! if compn wrong, 0 else 
                    ! its media index is obtained
      call epchkcmpn('epqCn2MediaIdx', compn)
      mediaidx = Det%Cn2media(compn)

      end

!        *************** inquire incident ptcl
      subroutine epqinc(aTrack)
      implicit none 
#include  "ZepTrackv.h"
       type(epTrack)::  aTrack
      aTrack =  Incident
      end
!     **************  inquire the first collision point
      subroutine epqFirstI(firstpos)
      implicit none
#include "ZepTrackv.h"
       type(epPos)::  firstpos  ! output.
      firstpos = FirstInt
      end
!     **************  inquire the first collision  track info.
      subroutine epq1stIntTrack(firstTrack)
      implicit none
#include "ZepTrackv.h"
       type(epTrack)::  firstTrack ! output.
      firstTrack = FirstIntTrack
      end

!     **************  inquire the comp. # where the first collision  (decay)
!                  etc occurred.   
      subroutine epqFirstICN(compno)
      implicit none
#include "ZepTrackv.h"
      integer,intent(out):: compno
      compno = FirstCn
      end
!     **************  inquire the first collision media
      subroutine epqFirstM(fmedia)
      implicit none
#include "Zcode.h"
#include "ZepTrackv.h"
       type(epmedia)::   fmedia  ! output
      if( Incident%p%code <= kmuon ) then
         if( FirstMedia%noOfElem == 1) then

!            for e/g/mu, collision element  is not fixed
!             but if the # of elements is only 1 we assign it
            FirstMedia%colElem = 1
            FirstMedia%colA = FirstMedia%A
            FirstMedia%colZ = FirstMedia%Z
         else
!           for e/g/mu  at present keep  colElem = 0 and colA,colA undef.
            FirstMedia%colElem = 0   ! for safety
         endif
      endif

      fmedia = FirstMedia

      end
      subroutine epqFirstP(proc)
      implicit none
#include "ZepTrackv.h"
      character*8  proc
      proc = Proc1
      end
      subroutine epqAlldE(alleloss)
!         returns sum of all energy loss in the detector
!       for one event (at end of 1event generation)
      implicit none
#include  "ZepTrackv.h"
      real(8),intent(out):: alleloss

      alleloss = sumdE
      end
