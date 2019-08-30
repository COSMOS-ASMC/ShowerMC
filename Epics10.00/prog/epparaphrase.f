      function epparaphrase(struc) result(ans)
      implicit none
#include "ZepMaxdef.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
!        This reurns detector structure paraphrased for 'new-1' etc.
!        paraphrasing is bi-dictional.  That is, if 'new-1' is 
!        paraphrased to 'cap', 'cap' can be paraphrased to 'new-1'
!        If no string matches, input is returned.
!
!     Before this is used,  epparadepo must be called to
!     store string pair to be praphrased.
! 
      character*(*)  struc  ! probably a max of MAX_STRUCCHR character string
!                           for which paraphrased string is to be found
!                           

      character(len=MAX_STRUCCHR) newpara(2, MaxNewStruc)

      character(*):: ans

      integer imax
      common /Zepparaphc/ newpara
      common /Zepparaph/ imax


 

      integer i 

      
      if(struc(1:4) .eq. 'new-') then
         do i = 1, imax
            if(newpara(1, i) .eq. struc) then
               ans  = newpara(2, i)
               return    ! **********
            endif
         enddo
         ans = trim(struc)
      else
         do i = 1, imax
            if(newpara(2, i) .eq. struc) then
               ans = newpara(1, i)
               return  ! ************
            endif
         enddo
         ans = trim(struc)
      endif
      end
      subroutine epiniparaph
      implicit none
      integer imax
      common /Zepparaph/ imax
      imax = 0
      end
!     ***************************
      function epparaphi(j) result(ans)
      implicit none
#include "ZepMaxdef.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
      integer j

      character(MAX_STRUCCHR):: newpara(2, MaxNewStruc)
      character(*):: ans
      integer imax
      common /Zepparaphc/ newpara
      common /Zepparaph/ imax


!        get j-th paraphrased string
      ans = newpara(2, j)
      end
!     ***************************        
      subroutine epparadepo(new, para)
      implicit none
#include "ZepMaxdef.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
      character*(*) new, para
!        deposit pair of string to be paraphrased.
      character(MAX_STRUCCHR):: newpara(2, MaxNewStruc)
      integer imax
      common /Zepparaphc/ newpara
      common /Zepparaph/ imax

      integer i

      if(imax .eq. 0) then
         do i = 1, MaxNewStruc
            newpara(1, i ) = ' '
            newpara(2, i ) = ' '
         enddo
      endif

      if(imax .lt. MaxNewStruc) then
         imax = imax +1
         if(new(1:4) .eq. 'new-') then
            newpara(1, imax) = new
            newpara(2, imax) = para
         elseif(para(1:4) .eq. 'new-') then
            newpara(1, imax) = para
            newpara(2, imax) = new
         else
            call cerrorMsg(new, 1)
            call cerrorMsg(para,1)
            call cerrorMsg(
     *     ' above two are  invalid for paraphrasing',0)
         endif
      else
         call cerrorMsg('too many paraphrasing data',0)
      endif
      end

!     ***********************************
      subroutine epparaphtbl(io)
      implicit none
#include "ZepMaxdef.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
      integer io ! input. logical dev. number of print
      character(MAX_STRUCCHR):: newpara(2, MaxNewStruc)
      integer imax
      common /Zepparaphc/ newpara
      common /Zepparaph/ imax
      integer i
!   
!           print paraphrase table
!

      if(imax .gt. 0) then
         write(io, *) ' ----------paraphrasing table---------'
      endif
      do i = 1, imax
         if(newpara(1, i) .ne. ' ') then
            write(io, *)
     *      i,'   ', newpara(1,i),' <--->  ',newpara(2,i)
         endif
      enddo
      end
!     *************************
      subroutine epqparatbl(tbl, maxp)
      implicit none
#include "ZepMaxdef.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
      character(MAX_STRUCCHR):: tbl(2, MaxNewStruc)
      integer maxp

      character(MAX_STRUCCHR):: newpara(2, MaxNewStruc)
      integer imax
      common /Zepparaphc/ newpara
      common /Zepparaph/ imax
!        returns newpara tble
      integer i
 
      maxp = imax
      do i = 1, imax
         tbl(1,i)= newpara(1,i)
         tbl(2,i)= newpara(2,i)
      enddo
      end
