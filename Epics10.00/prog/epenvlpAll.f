#include "ZsaveStruc.h"
!           ****************************************************
!
!         set the form parameter by seeing the configuration
!    form='box', 
!         xxxxx 'cyl','pipe', 'mix' will be set. xxxx
!         'box' may be set only if all components are box 
!               and all edges are parallel and have the  same
!               size in x, y. 
!     For other cases, 'mix' is set.
!         
!     If form != 'mix', faster execution will be possible.
!
        subroutine epsetform
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!     
        integer i

        do i = 1, Det%nct
           if(Det%cmp(i)%NMatreska .gt. 0) then
              form = 'mix'
              goto 10
           endif
        enddo
        call epcbox
 10     continue
        end
        subroutine epcbox
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include  "ZepPos.h"

        character*(*) formx

        character*(*) wstruc

        integer nw
       type(epPos)::  org, attr
!
        integer i
        real*8 ref

        
        ref = Det%cmp(1)%orgx
        do i = 1, Det%nct
           if(Det%cmp(i)%struc(1:3) .ne. 'box') goto 10
           if(Det%cmp(i)%orgx .ne. ref) goto 10
        enddo
        ref = Det%cmp(1)%orgy
        do i = 2, Det%nct
           if(Det%cmp(i)%orgy .ne. ref) goto 10
        enddo
        ref = Volat( Det%cmp(1)%vol+boxa)
        do i = 2, Det%nct
           if(Volat( Det%cmp(i)%vol+boxa ) .ne. ref) goto 10
        enddo
        ref = Volat( Det%cmp(1)%vol+boxb )
        do i = 2, Det%nct
           if(Volat( Det%cmp(i)%vol+boxb ) .ne. ref) goto 10
        enddo
        form = 'box'
        return         ! *************
 10     continue
        form = 'mix'
        return
!       *****************
        entry epqfrm(formx)
!           inquire form
            formx=form
            return

!       *****************
        entry  epqworld(nw, wstruc)
        nw = Det%nworld
        wstruc = Det%cmp(Det%nct)%struc
        return

!       ***********************
        entry epqwcoord(org, attr)
!          returns world origin and attribute
        org%x = Det%cmp(Det%nct)%orgx
        org%y = Det%cmp(Det%nct)%orgy
        org%z = Det%cmp(Det%nct)%orgz
        if(Det%cmp(Det%nct)%struc .eq.'box_w') then
!(((((((((((
           attr%x =Volat(Det%cmp(Det%nct)%vol+boxa )
           attr%y =Volat(Det%cmp(Det%nct)%vol+boxb )
           attr%z =Volat(Det%cmp(Det%nct)%vol+boxc )
!))))))))
        elseif( Det%cmp(Det%nct)%struc == 'sphere_w') then
!((((((
           attr%x = Volat( Det%cmp(Det%nct)%vol+sphr )
!))))))
        elseif( Det%cmp(Det%nct)%struc == 'cyl_w' .or.
     *    Det%cmp(Det%nct)%struc == 'cyl_y_w' .or.
     *    Det%cmp(Det%nct)%struc == 'cyl_x_w' .or.
     *    Det%cmp(Det%nct)%struc == 'cyl_z_w' ) then
           attr%x = Volat(Det%cmp(Det%nct)%vol+cylr)
           attr%y = Volat(Det%cmp(Det%nct)%vol+cylh)
        else
           write(0,*) ' world struc=', Det%cmp(Det%nct)%struc 
           write(0,*) ' strange '
           stop
        endif
        end
!     *************
      subroutine epqenvlper(i, org, abc)
!       inqure the enveloper of a given comp.
!       in world  coordinate
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

      integer i  ! input. component number 

       type(epPos)::   org  ! output. world orgin of a box. which
                           !    envelops  the i-th comp.
                           !  
       type(ep3Vec)::  abc  ! output. box a, b, c

      call epqenvlper0(i, org, abc)
      call epqenvlper1(i, org, abc)
      end      subroutine epqenvlper
!     *************
      subroutine epqenvlper0(i, org, abc)
!       inqure the enveloper of a given comp.
!       local coord.
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

      integer i  ! input. component number 

       type(epPos)::   org  ! output. local  orign of a box. which
                           !    envelops  the i-th comp.
                           !  
       type(ep3Vec)::  abc  ! output. box a, b, c


      character*80 msg
      character(len=MAX_STRUCCHR)::  epparaphrase, tempph
      character(len=MAX_STRUCCHR)::  basename


      integer::uscl

      if(i .le. 0 .or. i .gt.  Det%nct) then
         write(msg, *) ' inquired comp. #=', i, 
     *     ' not exists: epqenvlper'
         call cerrorMsg(msg, 0)
      endif
      call epGetBaseStrucName(Det%cmp(i)%struc, basename)
      if(basename .eq. 'box') then
         call epenvlpBox(Det%cmp(i), org, abc)

      elseif(basename .eq. 'cyl') then
         call epenvlpCyl(Det%cmp(i), org, abc)

      elseif(basename .eq. 'pipe') then
         call epenvlpPipe(Det%cmp(i), org, abc)

      elseif( basename .eq. 'prism' ) then
         call epenvlpPrism(Det%cmp(i), org, abc)
      elseif( basename  .eq. 'sphere') then
         call epenvlpSphere(Det%cmp(i), org, abc)
      else
         call epseeUnderScore(Det%cmp(i)%struc, uscl)
         tempph = epparaphrase(Det%cmp(i)%struc(1:uscl))
         if(tempph(1:4) .eq. 'new-') then
            call epenvlpNew(Det%cmp(i), org, abc)
         else
            write(msg,*)
     *      'struc=', Det%cmp(i)%struc,
     *      ' not supported: epqenvlper0. '
            call cerrorMsg(msg,0)
         endif
      endif

      end subroutine epqenvlper0

      subroutine epqenvlper1(i, org, abc)
!       convert output from epqenvloper0. 
!       in world  coordinate
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"

      integer i  ! input. component number 

       type(epPos)::   org  !in/output.  output from  eqpenvlper0.
                           !  world orgin of a box. which
                           !  envelops  the i-th comp.
                           !  
       type(ep3Vec)::  abc  ! in/output.output from  eqpenvlper0.
                           !   box a, b, c

       type(epPos)::  posw, maxp, orgt

      integer j, jmx
      real*8 dx(8), dy(8), dz(8)
      data dx/0., 1., 0., 1., 0., 1., 0., 1./
      data dy/0., 0., 1., 1., 0., 0., 1., 1./
      data dz/0., 0., 0., 0., 1., 1., 1., 1./
      if(NVTX .eq. 0) then
         jmx = 8
      else
         jmx = NVTX
      endif
      do j = 1, jmx
         if(NVTX .eq. 0) then
!             for all 8 corners of the box
!            posw.x = org.x + dx(j) * abc.x
!            posw.y = org.y + dy(j) * abc.y
!            posw.z = org.z + dz(j) * abc.z
            posw = epPos( org%x + dx(j) * abc%x,
     *                    org%y + dy(j) * abc%y,
     *                    org%z + dz(j) * abc%z)
         else
!            posw.x = org.x + VTXx(j) 
!            posw.y = org.y + VTXy(j) 
!            posw.z = org.z + VTXz(j)
            posw = epPos( org%x + VTXx(j), 
     *                    org%y + VTXy(j), 
     *                    org%z + VTXz(j) )
         endif
!            each vetex must be converted to world
!          to get the largest box
         call epl2w(i, posw, posw)
         if(j .eq. 1) then
            orgt = posw
            maxp = posw
         else
            orgt%x = min(orgt%x, posw%x)
            orgt%y = min(orgt%y, posw%y)
            orgt%z = min(orgt%z, posw%z)

            maxp%x = max(maxp%x, posw%x)
            maxp%y = max(maxp%y, posw%y)
            maxp%z = max(maxp%z, posw%z)
         endif
      enddo
      org = orgt
      abc%x = maxp%x -  org%x
      abc%y = maxp%y -  org%y
      abc%z = maxp%z -  org%z
      end  subroutine epqenvlper1

!     *************
      subroutine epenvlpAll
!          compute over all boundary of the system
!         in world  coordinate
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"


       type(epPos)::   orgx
       type(ep3Vec)::  abcx

       type(epPos)::  Orgsave, org
       type(ep3Vec):: Abcsave, abc
#ifdef USESAVE
      save   Orgsave, Abcsave
#endif
      integer i, imax

      Orgsave%x = 1.d10
      Orgsave%y = 1.d10
      Orgsave%z = 1.d10
      Abcsave%x = -Orgsave%x
      Abcsave%y = -Orgsave%y
      Abcsave%z = -Orgsave%z

      if( Det%nworld  .gt.  0 ) then
!          If the world has some attribute, it should be recognized
         if( Volat( Det%cmp(Det%nct)%vol + 1) .eq. 0.) then
!              no attribute so ignore world size
            imax = Det%nct -1
         else
!              see world size
            imax = Det%nct
         endif
      else
         imax = Det%nct
      endif
!((((((((((
!c      do   i=1, Det.nct
!))))))))
      do   i=1, imax
!           ignore component which is partially contained by another volume.
         if(Det%cmp(i)%NPContainer .eq. 0) then
!              i-th comp is not a partially contained one
            call epqenvlper(i, org, abc)
            call epLargerEnv(org, abc, Orgsave, Abcsave)
!                       here Abcsave is not yet edge length
!                    simeply max pos.
         endif
      enddo
!cccc  2002.02.09
!         Abc next is edge length 
      Abcsave%x = Abcsave%x - Orgsave%x
      Abcsave%y = Abcsave%y - Orgsave%y
      Abcsave%z = Abcsave%z - Orgsave%z
!cccc
      return
!     ************  inquire configuration
      entry epqcnf(orgx, abcx)
!     ******************
      orgx = Orgsave
      abcx = Abcsave
      end
!     **********************
      subroutine epLargerEnv(org, abc, orgx, maxpos)
!      suppose a virtual box's world origin be
!      'orgx' and it's largest (x,y,z) coordinate
!      is 'maxpos'. (maxpos = orgx + 3edge length)
!
!    This subroutine compares the box specified by org and abc
!    and returns a box which envelops both of them.
!    The new box is given by orgx and maxpos.
!
!
!       
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"


       type(epPos)::   org     ! input. world origin of a box

       type(ep3Vec)::  abc     ! input. abc of the box

       type(epPos)::  orgx     ! input/output ; in world
                              !  coord. 
       type(epPos)::  maxpos   ! input/output.; in world
                              !  coord.

       type(epPos)::   posw

      integer i, imx
      real*8 dx(8), dy(8), dz(8)
      data dx/0., 1., 0., 1., 0., 1., 0., 1./
      data dy/0., 0., 1., 1., 0., 0., 1., 1./
      data dz/0., 0., 0., 0., 1., 1., 1., 1./

      if(NVTX .eq. 0) then
         imx = 8
      else
         imx = NVTX
      endif
      do i = 1, imx
         if(NVTX .eq. 0) then
!             for all 8 corners of the box
            posw%x = org%x + dx(i) * abc%x
            posw%y = org%y + dy(i) * abc%y
            posw%z = org%z + dz(i) * abc%z
         else
            posw%x = org%x + VTXx(i) 
            posw%y = org%y + VTXy(i) 
            posw%z = org%z + VTXz(i)
         endif
!             now posw is one of imx corners in world coord.
         orgx%x = min(orgx%x, posw%x)
         orgx%y = min(orgx%y, posw%y)
         orgx%z = min(orgx%z, posw%z)
!             get max pos in world coord.
         maxpos%x = max(maxpos%x, posw%x)
         maxpos%y = max(maxpos%y, posw%y)
         maxpos%z = max(maxpos%z, posw%z)
      enddo
            
      end
!     *******************************
      subroutine epenvlpBox(comp, org, abc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"


       type(Component)::  comp
       type(epPos)::  org
       type(ep3Vec)::  abc

      org%x = 0.
      org%y = 0.
      org%z = 0.
!((((((((((
      abc%x = Volat( comp%vol+boxa )
      abc%y = Volat( comp%vol+boxb )
      abc%z = Volat( comp%vol+boxc )
!))))))))))
      NVTX = 0
      end

!     *******************************
      subroutine epenvlpCyl(comp, org, abc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"


       type(Component)::  comp
       type(epPos)::  org
       type(ep3Vec)::  abc
      real(8):: temp
!((((((((((((
      org%x = -Volat( comp%vol+cylr )
      org%y = org%x
      org%z = 0.
      abc%x = Volat( comp%vol+ cylr )*2
      abc%y = abc%x
      abc%z = Volat( comp%vol+cylh )
!))))))))))
      NVTX = 0
      if(comp%struc(1:5) == "cyl_y") then
         call epc2v_cyl(comp, org, org)
!         temp = abc.y
!         abc = ep3Vec(abc.x, abc.z, temp)
!         abc = ep3Vec(temp, abc.z, abc.x)
         call epc2v_cyl(comp, abc, abc)
      elseif(comp%struc(1:5) == "cyl_x") then
         call epc2v_cyl(comp, org, org)
!         temp=abc.x
!         abc = ep3Vec(abc.z, abc.y, temp)
!         abc = ep3Vec(abc.z,  temp, abc.y)
         call epc2v_cyl(comp, abc, abc)
      endif

      end

!     *******************************
      subroutine epenvlpPipe(comp, org, abc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"



       type(Component)::  comp
       type(epPos)::  org
       type(ep3Vec)::  abc

      real(8):: temp
!((((((((((
      org%x = -Volat( comp%vol+pipeor )
      org%y = org%x
      org%z = 0.
      abc%x = Volat( comp%vol+pipeor )*2
      abc%y = abc%x
      abc%z = Volat( comp%vol+ pipeh )
!))))))))))
      NVTX = 0


      if(comp%struc(1:6) == "pipe_y") then
         call epc2v_pipe(comp, org, org)
!         temp = abc.y
!         abc = ep3Vec(abc.x, abc.z, temp)
!         abc = ep3Vec(temp, abc.z, abc.x)
         call epc2v_pipe(comp, abc, abc)
      elseif(comp%struc(1:6) == "pipe_x") then
         call epc2v_pipe(comp, org, org)
!         temp=abc.x
!         abc = ep3Vec(abc.z, abc.y, temp)
!         abc = ep3Vec(abc.z, temp, abc.y)
         call epc2v_pipe(comp, abc, abc)
      endif
      end


!     *******************************
      subroutine epenvlpSphere(comp, org, abc)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepPos.h"


       type(Component)::  comp
       type(epPos)::  org
       type(ep3Vec)::  abc
!(((((((((
      org%x = -Volat( comp%vol+sphr )
!)))))))))
      org%y = org%x
      org%z = org%x

      abc%x = Volat( comp%vol+sphr ) *2
      abc%y = abc%x
      abc%z = abc%x
      NVTX = 0
      end


