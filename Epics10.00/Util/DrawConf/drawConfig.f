#include "ZepicsBD.h"
      program drawConfig
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include  "ZepDraw.h"


      integer menu, icon
      call copenf(11, '/tmp/$USER/Work/.configrc', icon)
      if(icon  .ne. 0) then
         write(0,*) ' config file not seen'
         stop
      else
         read(11, '(a)') configFile
         close( 11 )
      endif

      call epRdNewConfFile
      call cerrorMsg(' ', 1)
      call cerrorMsg(' ', 1)
      call cerrorMsg(' SELECT MENU NUMBER ', 1)

      do while (.true.)
         call cerrorMsg(' ', 1)
         call epMenu(menu)
         if(menu .eq. 0) then
            stop
         endif
         do while(.true.)
            if(menu .eq. 1) then
               if(subdnumber .eq. 0) then
                  call epExecGene(Det%cmp, Det%nct)
               else
                  call epExecGene(SubdArea(SubD(subdnumber)%loc+1), 
     *                            SubD(subdnumber)%nct)
               endif
            elseif(menu .gt. 1 .and. menu .le. maxmenu) then
               if(subdnumber .eq. 0) then
                  call epsetCond(Det%cmp,  menu)
               else
!                  call epsetCond(SubD(subdnumber).cmp, menu)
                  call epsetCond(
     *              SubdArea( SubD(subdnumber)%loc+1 ), menu)
               endif
            endif
            call cerrorMsg('Enter return (or take thin menu)', 1) 
            read(*,'(i2)', err=10) menu
            if(menu .eq. 0) goto 10
         enddo
 10      continue
      enddo
      end
!     **************
      subroutine epRdNewConfFile
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
      integer sum1, i, sum2,  levelmx

      call epfixConfig
      call epMemoCnfgFile   ! next one is here bef. 9.15
      call epReadConfig
      sellist = ' '
      call epsetDflt     ! set default parameters
      call epsetDflt2
      sum1 = 0
      do i = clevelmin, clevelmax
         sum1 = levelc(i) + sum1
      enddo

      sum2 = 0
      do i = 0, maxmaxlevel
         if(levelc(i) .gt. 0) levelmx=i
         sum2 = levelc(i) + sum2
      enddo

      if(sum1 .gt. maxdisplay) then
         call cerrorMsg(' ', 1)
         call cerrorMsg('**************************',1)
         write(msg,*)
     *   ' The  number of of components for display = ',
     *   sum1, ' is  large'
         call cerrorMsg(msg, 1)
         write(msg, *)
     *   ' For a slow computer, this may be a burden'
         call cerrorMsg(msg, 1)
         write(msg, *)
     *   ' If so, you may restrict the number of components for display'
         call cerrorMsg(msg, 1)
         call cerrorMsg(' by the levels of the comps.', 1)
            call cerrorMsg(' or by other restrictions ', 1)
         call cerrorMsg('****************************',1)
         call cerrorMsg('Enter return',1)
         read(*,*)
      elseif(sum2 .gt. maxdisplay .or. 
     *         levelmx .gt. clevelmax) then
         call cerrorMsg(' ', 1)
         call cerrorMsg(
     *   '****************************************************',1)
         call cerrorMsg(
     *  'You have a large number of components for display', 1)
         call cerrorMsg(
     *  'It may be a burden for a slow computer, so by default',1)
         call cerrorMsg(
     *    'some of the components with higher levles will not ',1)
         call cerrorMsg('be displayed.', 1)
         call cerrorMsg(
     *   'You can control the levels for display by menu 14',1)
         call cerrorMsg(
     *   '****************************************************',1)
         call cerrorMsg(' ', 1)
         call cerrorMsg('Enter return', 1)
         read(*,*)
      endif

      end

!     ***********************
      subroutine epMenu(menu)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
   
      integer menu  ! output.


      if(subdnumber .gt. 0) then
         call cerrorMsg('      ** SubDetector mode ** ', 1)
         write(msg, *) ' Selected sub-detector number=', subdnumber,
     *     ':  ', SubDName(subdnumber)
         call cerrorMsg(msg, 1)
      endif

      call cerrorMsg(' 1) Make data for display (=1)', 1)
      call cerrorMsg(' 2) Specify transparent surface/ranges', 1)
      call cerrorMsg(' 3) Specify components by number',1)
      call cerrorMsg(' 4) Specify components by media',1)
      call cerrorMsg(
     *     ' 5) Specify components by mother/daughter relation ', 1)
      call cerrorMsg(' 6) Specify components by structure ', 1)
      call cerrorMsg(' 7) Reread current config file', 1)
      call cerrorMsg(' 8) Show config info', 1)
      call cerrorMsg(' 9) Resume default setting',1)
      call cerrorMsg('10) Read new config data file',1)
      call cerrorMsg('11) Media<->color mapping', 1)
      call cerrorMsg(
     * '12) Save/recall display condition/change default',1)
      call cerrorMsg('13) Reverse Z-coord', 1)
      call cerrorMsg('14) Level control ', 1)
      if(NsubD .gt. 0) then   
         call cerrorMsg('15) Enter/Leave subdetector mode ', 1)
      endif
      call cerrorMsg(' 0) Exit', 1)
      menu = 1
      read(*,*)  menu
      end
!     ************************
      subroutine epsetCond(cmp, menu)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

       type(Component)::  cmp(*)
      integer menu

      if(menu .eq. 2) then
         call epRange
      elseif(menu .eq. 3) then
         call epbyNum
      elseif(menu .eq. 4) then
         call epbyMedia
      elseif(menu .eq. 5) then
         call epbyContain(cmp)
      elseif(menu .eq. 6) then
         call epbyStruc
      elseif(menu .eq. 7) then
         call epCloseDraw
         call epReadConfig
         if(subdnumber .gt. 0) then
            call epsetSubD2
         endif
      elseif(menu .eq. 8) then
         call epaboutFile
      elseif(menu .eq. 9) then
         call epsetDflt
         call epsetDflt2
      elseif(menu .eq. 10) then
         call epCloseDraw
         call epRdNewConfFile
      elseif(menu .eq. 11) then
         call epsetColMap
      elseif(menu .eq. 12) then
         call epdispCond
      elseif(menu .eq. 13) then
         call epreverse
      elseif(menu .eq. 14) then
         call eplevel
      elseif(menu .eq. 15) then
         call epsetSubD
      endif
      end
      subroutine epaboutFile
      implicit none
#include "ZepDraw.h"
       msg = 'Current config file: '//configFile
       call cerrorMsg(msg, 1)
       call cerrorMsg('Enter return',1)
       read(*,*)
       call epOutCnf(int(clevelmin), int(clevelmax), 6)
       end
!     ******************
      subroutine epRange
      implicit none
#include "ZepDraw.h"
      character*6 digits
      character*2 dig2
      integer keephow

      if(thetamin .gt. 360.) thetamin = thetamin - 360.
      write(msg,'(a, f6.1,a,f6.1,a)' ) 
     * 'Transparent azimuthal angle range for round objects'//
     * '(', thetamin,',', thetamax,') deg'
      call cerrorMsg(msg, 1)
      read(*,*) thetamin, thetamax
      if(thetamin .le. thetamax) then
         thetamin = thetamin + 360.
      endif
      write(msg, '(a,f6.1,a,f6.1,a)')
     * 'Opaque polar angle range of a'//
     * ' spherer (', pamin,',', pamax, ') deg'
      call cerrorMsg(msg, 1)
      read(*,*) pamin, pamax
      write(digits, '(i6)') how
      call kreplcha(digits, ' ', '0')
      keephow = how
 10   continue
      write(msg, *)
     *  'Enter 6 digits to specify surfaces to be '//
     *  'drawn for each box(=',digits,')'
      call cerrorMsg(msg, 1)
      call cerrorMsg('---put -1 if you need explanation---',1)
      read(*,*) how
      if(how .lt. 0) then
         call epboxsurf
         how = keephow
         goto 10
      endif
!         for cyl. like object
      write(dig2, '(i2)') howcyl
      call kreplcha(dig2, ' ', '0')
 20   continue
      write(msg, *)
     *  'Similarily enter 2 digits for floor and ceil of'//
     *  '  cylinder-like object(=',dig2,')'
      call cerrorMsg(msg, 1)
      read(*,*) howcyl

      end
!     ***************************
      subroutine epboxsurf
      implicit none
!       give explanation of box surface.
      call cerrorMsg(
     * "Suppose a canonical box with edge lengths 'a,b,c'", 1)
      call cerrorMsg('surface 1 is x-y at z=0', 1)
      call cerrorMsg('        6 is x-y at z=c', 1)
      call cerrorMsg('        2 is x-z at y=0', 1)
      call cerrorMsg('        5 is x-z at y=b', 1)
      call cerrorMsg('        3 is y-z at x=0', 1)
      call cerrorMsg('        4 is y-z at x=a', 1)
      call cerrorMsg('For example, 6 digits, 010001, specifies', 1)
      call cerrorMsg('surface 1 and 5 be drawn', 1)
      call cerrorMsg('Enter return',1)
      read(*,*)
      end
!     *******************************
      subroutine kreplcha(text, cha, by)
      implicit none
      character*(*) text  ! in/out.
      character*1 cha     ! in.
      character*1 by      ! in. All 'cha' in 'text' is replaced by 'by'

      integer i
      i = 1
      do while(i .ne. 0)
         i = index(text, cha)
         if(i .gt. 0) then
            text(i:i) = by
         endif
      enddo
      end
!     ****************************
      subroutine epbyNum
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

      integer i


      drawmode = 0
      do while(drawmode .ne. 1 .and. drawmode .ne. 2)
         write(msg,*) '   Total number of comps=', NumComp
         call cerrorMsg(msg, 1)
         call cerrorMsg(' 1) Draw only some components',1)
         call cerrorMsg(' 2) Hide only some components',1)
         call cerrorMsg(' Enter 1 or 2 (=1)',  1)
         drawmode = 1
         read(*,*) drawmode
      enddo
      call cerrorMsg(
     *   'Enter component numbers to be drawn or hidden with /', 1)
      call cerrorMsg(
     *   '1 2 5 -20 25 / means 1 2 5 6...19 20 25 /', 1)
      do i = 1,  maxspecifiable
         compnumb(i) = 0
      enddo
      read(*,*) compnumb
      call epsetnumcon(compnumb, drawmode, byNumb)  ! set drawor
      if(.not. byNumb) then
         call cerrorMsg('No comp. found for display',1)
         call cerrorMsg('condition cancelled', 1)
         do i = 1, NumComp
            draworhide(i) = 1
         enddo
      endif
      end
!     *********************
      subroutine epsetnumcon(numlist, dmode, yes)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
!          draworhide && numlist -> drawmode
      integer*4 numlist(*)
      integer dmode  ! in.  1 or 2. 1--> incl 2--> excl.
      logical yes   ! output. if there is some comp. to be shown
                    !         yes =.true.
      integer i
      
      logical eplisted, temp

      yes = .false.
      do i = 1, NumComp
         temp = eplisted(i, numlist, maxspecifiable)
         if( (temp .and. dmode .eq. 1)  .or.
     *       (.not. temp .and. dmode .eq. 2 ) ) then
            draworhide(i) = 1 
            if(.not. yes) yes = draworhide(i) .gt. 0
         else
            draworhide(i) = 0
         endif
      enddo

      end
!     ***********************************
      logical function eplisted(i, numlist, n)
      implicit none
      integer i   ! input.  see if "i" is listed in numlist
      integer n   ! input.  size of numlist
      integer*4 numlist(n)  ! input. if 0 appears, no more data is assumed

      integer j
      do j = 1, n
         if(numlist(j) .eq. 0) goto 5
         if(numlist(j) .eq. i) then
            eplisted = .true.
            goto 10
         endif
         if(numlist(j) .lt. 0) then
            if(j .gt. 1) then
               if(abs(numlist(j-1)) .le. i .and.
     *           i .le. abs(numlist(j)) ) then
                  eplisted = .true.
                  goto 10
               endif
            endif
         endif
      enddo
 5    continue 
      eplisted = .false.
 10   continue
      end
!     ********************
      subroutine epbyMedia
      implicit none
#include "ZepDraw.h"
      
      integer klena, i, j

      call cerrorMsg('Available media is as follows', 1)
      call cerrorMsg(mlist, 1)
      if(sellist .eq. ' ' ) then
         sellist = mlist
      endif
      
      selmedia = 0
      do while (selmedia .le. 0 .or. selmedia .gt. 2)
         call cerrorMsg(
     *   ' 1) Select some media to be shown',1)
         call cerrorMsg(
     *   ' 2) Select some media to be hidden',1)
         call cerrorMsg('Enter 1 or 2', 1)
         read(*,*) selmedia
      enddo
      write(msg, *)
     * "Enter media name list to show or to hide"//
     * "(case sensitive; like: "// 
     *  mlist(1:klena(mlist))//")"
      call cerrorMsg(msg, 1)
      sellist = mlist
      read(*, '(a)') sellist

      call kgetField(sellist, seledMedia, maxsel, nselMedia)
      byMedia = .false.

      do i = 1, nselMedia
         do j = 1, noOfMedia
            if(seledMedia(i) .eq. matdef(j)) then
               byMedia = .true.
               exit
            endif
         enddo
      enddo
 10   continue
      if(.not. byMedia) then
         call cerrorMsg('Your input not exist in the mediia list.',1)
         call cerrorMsg('media selection cancelled',1)
      endif
      end
!     ****************
      subroutine epbyContain(cmp)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include "ZepDraw.h"

       type(Component)::  cmp(*)

      integer sel
      sel = -1
      do while( ( sel .lt. 0 .or. sel .gt. 3) )
         call cerrorMsg(
     *  ' 1) Draw mothers of specified daughters',1)
         call cerrorMsg(
     *  ' 2) Draw daughters of spcified mothers',1)
         call cerrorMsg(
     *  ' 3) show mother/daughter component list', 1)
         call cerrorMsg('Select number(=1)', 1)
         sel = 1
         read(*,*) sel
      enddo

      if(sel .ge. 1 .and. sel .le. 2) then
         containmode = sel
         call epsetContain(cmp)
      elseif(sel .eq. 3) then
         call epshowMother(cmp)
         call epshowDaughter(cmp)
      endif
      end
!     **********************
      subroutine epsetContain(cmp)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
      
       type(Component)::  cmp(*)

      integer i
      
!      call cerrorMsg(
!     *  'Note: daughters and mothers are both drawn. ', 1)
      call cerrorMsg(
     *  "A negative number means 'partially containing'", 1) 

      if(containmode .eq. 1) then
         call cerrorMsg(
     *   'Enter daughter comp. numbers to'//
     *     ' display their mothers (with / last )', 1)
      else
         call cerrorMsg(
     *   'Enter mother comp. numbers to'//
     *   ' display their daughters (with / last )', 1)
      endif

      do i = 1,  maxspecifiable
         mamordaught(i) = 0
      enddo
      read(*,*) mamordaught
      call epifanybyContain(cmp)
      if(.not. byContain) then
         call cerrorMsg(
     *  'No comp. to be drawn by your condition;'//
     *  ' The condition cancelled', 1)
      endif
      end
!     ***************************
      subroutine epifanybyContain(cmp)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

       type(Component)::  cmp(*)

      integer i, j, nm, mc, k
!          see if at least on component exists to display
!      if not, set byContain =.false.
      byContain = .false.
      if(containmode .eq. 1) then
!           search a mother of specified daughters
         do i = 1, maxspecifiable
            j = mamordaught(i) 
            if(j .eq. 0) goto 100
            do k = 1, NumComp
               nm =cmp(k)%NMatreska
               do mc = 1, nm
! ((((((((((((((((((
!                  if(cmp(k).ContainsR(mc) .eq. j)  then
                  if(CnArea( cmp(k)%ContainsR+mc ) .eq. j)  then
!)))))))))))0
                     byContain = .true.
                     goto 100
                  endif
               enddo
            enddo
         enddo
      else
!         search a daughter of a given mother
         do i = 1, maxspecifiable
            j = mamordaught(i)
            if(j .eq. 0) goto 100
!             |j| is a mother. j < 0 means that daugther is partially
!                  contained by mother, |j|.
            if(abs(j) .lt. NumComp) then ! must < 
               if(cmp(abs(j))%NMatreska .gt. 0) then
                  nm = cmp(abs(j))%NMatreska 
                  if(j .gt. 0) then
                     byContain = .true.
                  else
                     do mc = 1, nm
!((((((((((((
!                        if(cmp(abs(j)).ContainsR(mc) .lt. 0) then
                        if( CnArea( cmp(abs(j))%ContainsR+mc )
     *                       .lt. 0) then
!))))))))))))
                           byContain = .true.
                           goto 100
                        endif
                     enddo
                  endif
               endif
            endif
         enddo
      endif
 100  continue
      end
!     ***************************
      subroutine epexContain(cmp, i, icon)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

       type(Component)::  cmp(*)

      integer i    ! input.  comp. number
      integer icon  ! output. =0, if i-th comp. satisies cond. for display
                    !             else 1.

      integer k, nm, mc, j

      icon = 1

!
      if(containmode .eq. 1) then
!          see if i-th comp. is a mother  of one of mamordaught(*)
!          if so ,make icon  = 0
!
         nm = cmp(i)%NMatreska
         do mc = 1, nm
            do k = 1, maxspecifiable
               if(mamordaught(k) .ne. 0) then
!(((((((((((((((((
!                  if( cmp(i).ContainsR(mc) .eq. mamordaught(k)) then
                  if( CnArea( cmp(i)%ContainsR+mc ) 
     *               .eq. mamordaught(k)) then
!)))))))))))
                     icon = 0
                     goto 100
                  endif
               else
                  goto 10
               endif
            enddo
 10         continue
          enddo
       else
!           see if i-th component is a daughter of one of mamordaught(*)
          do k = 1, maxspecifiable
             j = mamordaught(k) 
             if(j .eq. 0) goto 100
             if(abs(j) .le. NumComp) then
                nm = cmp(abs(j))%NMatreska
                do mc = 1, nm
                   if( (j .gt. 0. .and.
!((((((((((((((((((((
!     *                cmp(abs(j)).ContainsR(mc) .eq. i) .or.
!     *                 ( j .lt. 0 .and.
!     *                 cmp(abs(j)).ContainsR(mc) .eq. -i) )  then
     *                CnArea( cmp(abs(j))%ContainsR+mc ) .eq. i) .or.
     *                 ( j .lt. 0 .and.
     *                CnArea(cmp(abs(j))%ContainsR+mc) .eq. -i) )  then
!)))))))))))))))
                      icon = 0
                      goto 100
                   endif
                enddo
             endif
          enddo
       endif
 100   continue
       end
             
!     **************************
      subroutine epaddifyet(a, n, j)
      implicit none
      integer a(*) ! in/out.
      integer n    ! in/out.  currently filled max 'a' pos.
      integer j    ! in.  if j  is not in a, j is added in a and n is
                   !          incremented.
      integer i
      do i = 1, n
         if(a(i) .eq. j) then
            goto 10
         endif
      enddo
      n = n + 1
      a(n) = j
 10   continue
      end

!     ********************
      subroutine epshowMother(cmp)
      implicit none
#include "ZepDraw.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"

       type(Component)::  cmp(*)

!               show Matreska component
      integer i, k, j, lc, ll, m

      do i = 1, (NumComp - NumWorld)
         k = cmp(i)%NMatreska 
         if( k .gt. 0) then
            ll = min(10, k)
            write(msg,
     *     '(" #",i7," contains " 10i7)') i,
     *      (CnArea( cmp(i)%ContainsR+j ),j=1, ll)
            call  ksupblank(msg, lc)
            call cerrorMsg(msg, 1)
            if(k > 10) then
               do m = 11, k, 10
                  write(msg, 
     *            '("     " 10i7)') 
     *            (CnArea( cmp(i)%ContainsR+j ),j=m, min((m+9), k))
                  call  ksupblank(msg, lc)
                  call cerrorMsg(msg, 1)
               enddo
            endif
         endif
      enddo
      end
!     ********************
      subroutine epshowDaughter(cmp)
      implicit none
#include "ZepDraw.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
       type(Component)::  cmp(*)
!            show components which are contained by others
      integer i, k, lc

      do i = 1, (NumComp - NumWorld)
         k = cmp(i)%NContainer
         if( k .gt. 0) then
            write(msg,
     *      '(" #",i7," is contained by " i7," component(s)")') i,
     *      k
            call  ksupblank(msg,  lc)
            call cerrorMsg(msg, 1)
         endif
      enddo
      end
!     ***********************
      subroutine epExecGene(cmp, ncmp)
!         epexec for subd
      implicit none
#include  "Zglobalc.h"
#include "ZepMaxdef.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include "ZepDraw.h"

      integer ncmp
       type(Component)::  cmp(ncmp)

      integer klena, i, icon
      logical first, showusage 
      character(len=MAX_STRUCCHR) epparaphrase, tempph
      integer  nv
      character*2  dql/' "'/, dqr/'" '/
      real*8 r
      data showusage/.false./  ! now alwasys.  actual usage will be shown
                    ! by script 
      save showusage
      integer uscl

!      call mediaCount(cmp, matdef, noOfMedia, compnToMn)
      call mediaCount(cmp, ncmp, matdef, noOfMedia)
      call epOpenMediaFile   ! open media file
      do i = 1, ncmp
         nv = 0
         call epaccept(cmp, i, icon)
         if(icon .ne. 0) goto 500
!           i is being drawn
         usedmedia(compnToMn(i)) = .true.
         if(cmp(i)%struc(1:3) .eq. 'box') then
            call epdrawBox(cmp(i), pv, nv)
         elseif(cmp(i)%struc(1:3) .eq. 'cyl') then
            call epdrawCyl(cmp(i), pv, nv)
         elseif(cmp(i)%struc(1:6) .eq. 'sphere') then
            r = Volat( cmp(i)%vol + sphr)
            call epdrawSphere(r,
     *        r*cos(pamax*Torad), r*cos(pamin*Torad), pv, nv)
         elseif(cmp(i)%struc(1:4) .eq. 'pipe') then
            call epdrawPipe(cmp(i), pv, nv)
         elseif(cmp(i)%struc .eq. 'prism' .or.
     *          cmp(i)%struc(1:6) .eq. 'prism_' ) then
            call epdrawPrism(cmp(i),  pv, nv)
         elseif(cmp(i)%struc(1:4) .eq. 'new-' ) then
            call epDrawNew(cmp(i), pv, nv)
         else
            call epseeUnderScore( cmp(i)%struc, uscl)
            tempph = epparaphrase(cmp(i)%struc(1:uscl))
            if(tempph(1:4) .eq. 'new-' ) then
               call epDrawNew(cmp(i), pv, nv)
            else
               call cerrorMsg(cmp(i)%struc(1:8), 1)
               call cerrorMsg(' is undefined', 0)
            endif
         endif
!              output into disk
         call epwvtx(i, compnToMn(i)+offset, cmp(i),
     *        nv,  pv)
         if(nv .gt. maxvtx) then
            call cerrorMsg('too many vertex in one comp.', 1)
            call cerrorMsg(cmp(i)%struc, 0)
         endif
 500     continue
      enddo
!     
      call copenfw(21, '/tmp/$USER/Work/gnu', icon)
      call copenfw(22, '/tmp/$USER/Work/gnu2',icon)
      write(21, '(a)') 'set para'
      write(22, '(a)') 'set para'
      write(21, '(a)') 'set nokey'
      write(22, '(a)') 'set nokey'
      write(21, '(a)') 'set hidden'
      write(22, '(a)') 'set nohidden'
     
      first = .true.
      do i =  1, noOfMedia
         if( usedmedia(i) ) then
            if(first) then
               write(21, '("splot ",a,a,a, " w l lc",i4)' )  
     *              dql, fn(i)(1:klena(fn(i) )), dqr, mcolor(i) 
               write(22, '("splot ",a,a,a, " w l lc ",i4)' )  
     *              dql, fn(i)(1:klena(fn(i) )), dqr, mcolor(i) 
               first = .false.
            else
               write(21, '("rep ",  a, a, a, " w l lc", i4)') 
     *              dql, fn(i)(1:klena(fn(i))), dqr, mcolor(i)
               write(22, '("rep ",a,a,a, " w l lc ",i4)' )  
     *              dql, fn(i)(1:klena(fn(i) )), dqr, mcolor(i) 
            endif
         endif
      enddo
!c      write(21, '(a)')'set nohidden'
      if(showusage) then
!                  ======= show usage ========
         call cerrorMsg(
     *   ' *********        Useage  **************', 1)
         call cerrorMsg(
     *   ' In gnuplot, you can draw the configuration by ', 1)
         call cerrorMsg(
     *   '   call "config"   (hidden mode; slow)', 1) 
         call cerrorMsg('OR',1)
         call cerrorMsg(
     *   '   call "xconfig"   (nohidden mode; faster)', 1) 
         call cerrorMsg(
     *   ' or  later by  "clear"; call "all" (etc);', 1)
         call cerrorMsg(
     *   ' (retry if all of the components not shown)', 1)
         call cerrorMsg(
     *   ' -----------------------------------------', 1)
         showusage = .false.
      else
         call cerrorMsg(
     *   'ready for config display by  geomview or gnuplot', 1)
      endif
      close(21)
      close(22)
      call epCloseDraw
      end
!     ******************
      subroutine epaccept(cmp, i, icon)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include "ZepDraw.h"

       type(Component)::  cmp(*)

      integer i     ! in:see if i-th comp. is to be displayed
      integer icon  ! out: if so, icon = 0, else 1



      if(draworhide(i) .eq.  0) then
         icon = 1
         goto 500  !  if some cond. not satified, no display
      else
         icon = 0
      endif
      if(byMedia) then
         call epexMedia(cmp, i, icon)
         if(icon .ne. 0) goto 500
      endif

      if(byStruc) then
        call  epexStruc(cmp, i, icon)

        if(icon .ne. 0) goto  500
      endif

      if(byContain) then
         call epexContain(cmp, i, icon)

         if(icon .ne. 0) goto 500
      endif

      if( cmp(i)%level .gt. clevelmax  .or.
     *    cmp(i)%level .lt. clevelmin ) then
         icon =1
         goto 500
      endif
 500  continue
      end
      subroutine epwvtx(ith, io, cmp,  n, p)
      implicit none
#include "ZepMaxdef.h"
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
#include "ZepDraw.h"

      integer ith   ! input. ith component number is drawn
       type(Component)::  cmp
      integer io            ! input.  logical device for output.
      integer  n            ! input.  number of vertex in p
       type(epPos)::   p(n)  ! input.  vertex for gnuplot splot
      character(len=MAX_STRUCCHR) type

       type(epPos)::  pw

      integer i


      do i = 1, n
         if(p(i)%x .eq. gpsep) then
            write(io, *) ' '
         else
!            call epl2w(cn, p(i), pw)  ! cannot be used because
                                       ! subdetectors may be drawn 
             call epl2wx(cmp, p(i), pw)  

             if(cmp%subdidx .ne. 0) then
                type = SubDName(cmp%subdidx)
             else
                type = cmp%struc
             endif

             if(Reverse .eq. 1 ) then
                write(io, '(3g13.5,i6," ", a, i5," ", a)')
     *          sngl(pw%x), sngl(pw%y), sngl(pw%z), 
     *          ith, cmp%struc, cmp%fsubdc, type
             else
                write(io, '(3g13.5, i6, " ", a, i5," ", a)' )
     *          sngl(pw%x), sngl(pw%y), sngl(Zrevmax-pw%z),
     *          ith, cmp%struc, cmp%fsubdc, type
             endif
         endif
      enddo
      if( n .gt.  0) then
         write(io, *)  '  '
      endif
      end
      subroutine epl2wx(cmp, pos, poso)
!         local to world conversion
!        slight modification of epl2w
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"

       type(Component)::  cmp
       type(epPos)::  pos ! input
       type(epPos)::  poso ! output


      if( cmp%rotation) then
         poso%x = cmp%direc(1) * pos%x +
     *      cmp%direc(4) * pos%y +
     *      cmp%direc(7) * pos%z


         poso%y = cmp%direc(2) * pos%x +
     *      cmp%direc(5) * pos%y +
     *      cmp%direc(8) * pos%z

         poso%z = cmp%direc(3) * pos%x +
     *      cmp%direc(6) * pos%y +
     *      cmp%direc(9) * pos%z

         poso%x = poso%x + cmp%orgx
         poso%y = poso%y + cmp%orgy
         poso%z = poso%z + cmp%orgz
      else
         poso%x = pos%x + cmp%orgx
         poso%y = pos%y + cmp%orgy
         poso%z = pos%z + cmp%orgz
      endif
      end
!     **************  
      subroutine epdispCond
      implicit none
#include "ZepDraw.h"
!
      integer which
      which = -1
      do while (which .le. -1 .or. which .ge. 6)
         call cerrorMsg('0--Cancel', 1)
         call cerrorMsg('1--Show current display condition',1) 
         call cerrorMsg(
     *   "2--Save current condition but don't change default",
     *    1)
         call cerrorMsg(
     *   "3--Save current condition and make it default",
     *    1)
         if(Nsaved .gt. 1) then
            call cerrorMsg('4--Recall saved condition', 1)
            call cerrorMsg('5--Change default conditon',1)
         endif
         read(*,*) which
      enddo
      if(which .eq. 1) then
         call epshowCond
      elseif(which .eq. 2) then
         call epsaveCond(0)
      elseif(which .eq. 3) then
         call epsaveCond(1)
      elseif(which .eq. 4) then
         call eprecallCond
      elseif(which .eq. 5) then
         call epchgDfltCond
      endif
      end
      subroutine epchgDfltCond
!        change default cond
      implicit none
#include "ZepDraw.h"
      if(Nsaved .gt. 1) then
         call epshowCond
         call cerrorMsg('choose default display cond.',1)
         read(*,*) Ndflt
         if(Ndflt .ge. 1 .and. Ndflt .le. Nsaved) then
            call epselCond(Ndflt)
            call eprwCond    ! rewrite  cond.
         endif
      endif
      end
!
      subroutine epsaveCond(info)
!        save status now
      implicit none
#include "ZepDraw.h"
      integer info  ! input.  =0 don't change defalut
                    !         =1 chande defalut to the new one
!         read saved cond
!       call epgetCond
!         show  it
!       call epshowCond
       if(Nsaved .ge. NsavedMX) then
          call cerrorMsg(
     *         'no more room to save display conditions', 1)
          call cerrorMsg(
     *     'Some must be deleted', 1)
          call epdelCond
       endif
       if(Nsaved .lt. NsavedMX) then
          if(info .ne. 0) then
             Ndflt = Nsaved + 1
          endif
          call epaddCond
       endif
       end
      subroutine eprecallCond
      implicit none
#include "ZepDraw.h"

      integer sel
!         read saved cond
      call epgetCond
!         show  it
      call epshowCond
      call cerrorMsg('Select nubmer or enter 0 to cancel' ,1)
      read(*,*) sel
      if(sel .ge. 1 .and. sel .le. Nsaved) then
         call epselCond(sel)
      endif
      end
      subroutine epselCond(sel)
      implicit none
#include "ZepDraw.h"
      integer sel
      how = ahow(sel)
      howcyl = ahowcyl(sel)
      byContain = abyContain(sel)
      byMedia = abyMedia(sel)
      byStruc = abyStruc(sel)
      thetamin = athetamin(sel)
      thetamax = athetamax(sel)
      pamin = apamin(sel)
      pamax = apamax(sel) 
      drawall= adrawall(sel)
      end

      subroutine epgetCond
      implicit none
#include "ZepDraw.h"

!         read saved cond
      integer icon
      integer i
!     
      call copenf(10, CondFile, icon)
      if(icon .eq. 0 ) then
         read(10, *)  Nsaved, Ndflt
         do i = 1, Nsaved
            read(10, '(a)') CondID(i)
            read(10, *)   ahow(i), ahowcyl(i), abyContain(i), 
     *           abyMedia(i),  abyStruc(i), 
     *           athetamin(i), athetamax(i), apamin(i), apamax(i), 
     *           adrawall(i)
         enddo    
         close(10)
      else
         call cerrorMsg("Don't worry the above message",1)
         Nsaved = 0
         Ndflt = 0
      endif

      end
      subroutine epshowCond
      implicit none
#include "ZepDraw.h"

!        show cond

      integer i
      if(Nsaved .ge. 1) then
         write(msg,*)
     *    'currently available display conditions:  Defalut=',
     *     Ndflt
         call cerrorMsg(msg, 1) 
      else
         call cerrorMsg(
     *   'current display condition is system default:',1)
         write(msg, *) 'box=', how,' cyl=', howcyl,
     *   ' bycontain=', byContain, ' bymedia=',
     *      byMedia, ' bystruc=', byStruc, 
     *     ' thetamin=', thetamin,
     *     ' thetamax=', thetamax,
     *     ' azimuthmin=', pamin,
     *     ' azimuthmax=', pamax,
     *     ' draw all=', drawall
         call cerrorMsg(msg, 1)
      endif
      do i = 1, Nsaved
         write(msg, '(i5, a, a)') i,' ****** ', CondID(i)
         call cerrorMsg(msg, 1)
         write(msg, *) 'box=', ahow(i),' cyl=', ahowcyl(i),
     *   ' bycontain=', abyContain(i), ' bymedia=',
     *      abyMedia(i), ' bystruc=', abyStruc(i), 
     *     ' thetamin=', athetamin(i),
     *     ' thetamax=', athetamax(i),
     *     ' azimuthmin=', apamin(i),
     *     ' azimuthmax=', apamax(i), 
     *     ' draw all=', adrawall(i)
         call cerrorMsg(msg, 1)
      enddo
      end
      subroutine epdelCond
      implicit none
#include "ZepDraw.h"

!         rm some cond
      integer i, j, k

      call epshowCond
      call cerrorMsg(
     *    'Enter cond number to be removed',1)

      read(*,  *) i

      if(i .ge. 1 .and. i .le. Nsaved) then
         k = 0
         do j = 1, Nsaved
            if(j .ne. i) then
               k = k +1
               if(k .ne. j) then
                  CondID(k) = CondID(j)
                  ahow(k) =  ahow(j)
                  ahowcyl(k) = ahowcyl(j)
                  abyContain(k) = abyContain(j)
                  abyMedia(k) = abyMedia(j) 
                  abyStruc(k) = abyStruc(j)
                  athetamin(k) = athetamin(j)
                  athetamax(k) = athetamax(j)
                  apamin(k) = apamin(j)
                  apamax(k) = apamax(j) 
                  adrawall(k)= adrawall(j)
               endif
            endif
         enddo
         Nsaved = k
      else
         call cerrorMsg('no such number', 1)
      endif               
      end
      subroutine epaddCond
      implicit none
#include "ZepDraw.h"
      integer icon, i

      call copenfw(12, CondFile, icon)
      if(icon .ne. 0) then
         call cerrorMsg('display condition cannot be saved', 1)
      endif
      if(Nsaved .ge. NsavedMX) then
         call cerrorMsg(
     *    'display condition cannot be saved; too many', 1)
      else 
         Nsaved = Nsaved + 1
         call cerrorMsg(
     *  'Enter some comment for this condition(<24 chrs)', 1)
         read(*,'(a)') CondID(Nsaved)

         write(12, *) Nsaved, Ndflt 

         do i = 1, Nsaved-1
            write(12, *) CondID(i)
            write(12, *)   ahow(i), ahowcyl(i), abyContain(i), 
     *           abyMedia(i),  abyStruc(i), 
     *           athetamin(i), athetamax(i), apamin(i), apamax(i), 
     *           adrawall(i), ' /'
         enddo    
         write(12, *) CondID(Nsaved)
         write(12, *) how, howcyl, byContain,
     *     byMedia,  byStruc, thetamin, thetamax, pamin, pamax, 
     *     drawall, ' /'
         close(12)
      endif
      end
      subroutine eprwCond
      implicit none
#include "ZepDraw.h"
      integer icon, i

      call copenfw(12, CondFile, icon)
      if(icon .ne. 0) then
         call cerrorMsg('display condition cannot be saved', 1)
      else
         write(12, *) Nsaved, Ndflt 
         do i = 1, Nsaved
            write(12, *) CondID(i)
            write(12, *)   ahow(i), ahowcyl(i), abyContain(i), 
     *           abyMedia(i),  abyStruc(i), 
     *           athetamin(i), athetamax(i), apamin(i), apamax(i), 
     *           adrawall(i), ' /'
         enddo    
         close(12)
      endif
      end
!     **************    
      subroutine epsetDflt
      implicit none
#include "ZepDraw.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"

      integer i, klena, j


      logical first/.true./

      character*6 gplotcol(10)
      data gplotcol/'black','green', 'blue', 'red', 'purple',
     *    'cyan', 'yellow','dblack', 'orange','gray'/

      save first

      do i = 1, 10
         origcol(i) = gplotcol(i)
      enddo
      subdnumber = 0
      if(first) then
         CondFile = 'Draw_Cond'
!         how =011100
         how = 111111
         Reverse = -1  !  this is reverted in epreverse.
         call epreverse
!         howcyl = 00
         howcyl = 11
         byContain = .false.
         byMedia = .false.
         byStruc = .false.
         thetamin = 260. + 360.
         thetamax = 310
         pamin = 10.
         pamax = 170.
         drawall = .true.
         clevelmax = 6
         clevelmin = 0
         call epgetCond         ! get old conditions, if any 
         first = .false.
      endif
      
      if(Ndflt .ge.1 .and. Ndflt .le. NsavedMX) then
         call epselCond(Ndflt)
      endif
!        default is to draw  all.
      do i = 1,  maxcomp
         draworhide(i) = 1
      enddo
      write(mlist, '(40a)')
     *   (matdef(i)(1:klena(matdef(i)))//' ', i=1,noOfMedia)

!         space is not drawn by default.
      sellist = 'sp'
      selmedia = 2
      call kgetField(sellist, seledMedia, maxsel, nselMedia)
      byMedia = .true.

      do i = 1, noOfMedia
         mcolor(i) = i-1
         j =mod(i-1,9)+1
         collist(i) =origcol(j)
      enddo
      end

      subroutine epsetDflt2
      implicit none
#include "ZepMaxdef.h"
#include "ZepDraw.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"

      integer i, j, klena
!         make structure table
      nstruc = 0
      do i = 1, Det%nct
         do j = 1, nstruc
            if(strucarray(j) .eq. Det%cmp(i)%struc) then
               goto 10
            endif
         enddo
         nstruc = nstruc + 1
         if( nstruc > maxNewStruc ) then
            write(0,*) ' # of new structure > ', maxNewStruc
            write(0,*) ' Enlarge MAX_NEW_STRUC in ZepMaxdef%h'
            stop
         endif
         strucarray(nstruc) =Det%cmp(i)%struc
         write(0,*) nstruc,  strucarray(nstruc)
 10      continue
      enddo
      write(struclist, '(40a)')
     *   (strucarray(i)(1:klena(strucarray(i)))//' ', i=1,nstruc)

      NumComp = Det%nct
      NumWorld = Det%nworld
      end
!     ******************
      subroutine epsetColMap
      implicit none
#include "ZepDraw.h"
!          change mcolor;
      integer i, j, k
      integer list(maxmat*2)
      write(msg, '(10(i4," ",a))') (i-1, origcol(i), i=1,10)
      call cerrorMsg(
     * 'avialable colors (mapping is VERY GNOMISH)',1)
      call cerrorMsg(msg,1)
      call cerrorMsg('current color mapping is:', 1)
      write(msg, '("media #   media    color # color")')
      call cerrorMsg(msg, 1)
      do i = 1, noOfMedia
!         j =mod(mcolor(i),9)+1
         write(msg, '(i6,"   ", a, i8,"   ",a)')
     *        i, matdef(i), mcolor(i),  collist(i)
         call cerrorMsg(msg, 1)
      enddo
      call cerrorMsg(
     *  'spcify mapping by a list of pair of numbers with last / ',1)
      call cerrorMsg(
     *  'E%g., 2 5 6 1/',1)
      call cerrorMsg(
     * ' whcih means 2nd meida by color number 5; 6th  by 1',1)
      do i =1, noOfMedia*2-1, 2
         list(i) =  -1
         list(i+1) = -1
      enddo
      read(*,*) list
      do i = 1, noOfMedia*2 -1, 2
         if(list(i) .eq. -1 .or. list(i+1) .eq. -1 ) goto 10
         j = list(i)
         mcolor( j ) = list(i+1)
         k = mod(mcolor(j),9 )  + 1
         collist(j) = origcol(k)
      enddo
 10   continue
      end
!     ******************
      subroutine epsetSubD
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

!         select sub detector

      integer i

      call cerrorMsg('  #  Subdetector name', 1) 
      do i = 1, NsubD
         write(msg, *) i,'    ', SubDName(i)
         call cerrorMsg(msg, 1)
      enddo

      subdnumber = -1
      do while (subdnumber .lt. 0 .or. subdnumber .gt. NsubD)
         call cerrorMsg('Select subdetector # or enter 0',1)
         read(*, *) subdnumber
      enddo
      if(subdnumber .eq. 0) then
         call epsetDflt2
      else
         call epsetSubD2
      endif
      end
!     ******************
      subroutine epsetSubD2
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

!         select sub detector

      integer i, j, klena

      NumComp = SubD(subdnumber)%nct
      NumWorld = SubD(subdnumber)%nworld
      nstruc = 0 
      do i = 1, NumComp
         do j = 1, nstruc
            if(strucarray(j) .eq. 
     *          SubdArea(SubD(subdnumber)%loc+i)%struc) then
               goto 10
            endif
         enddo
         nstruc = nstruc + 1
         strucarray(nstruc) =
     *        SubdArea( SubD(subdnumber)%loc+i )%struc
 10      continue
      enddo
      write(struclist, '(40a)')
     *        (strucarray(i)(1:klena(strucarray(i)))//' ', i=1,nstruc)
      end
!    *******************
      subroutine epreverse
      implicit none
#include "ZepDraw.h"
       type(epPos)::  org, abc
      integer icon, lc

      if(Reverse .eq. 1) then
         Reverse = -1
         call epqcnf(org, abc)
         Zrevmax = org%z + abc%z
      else
         Reverse = 1
         Zrevmax = 0.
      endif
!          following is for drawing cascade trace z-reversed
      call copenfw(12, '.Reverse', icon)
      if(icon .ne. 0) then
         call cerrorMsg('.Reverse file cannot be opened', 0)
      endif
      write(msg, *) Reverse
!        revmove blank at head ( gnuplot awk problem )
      call kseblk(msg, ' ', lc)      
      write(12,'(a)') msg(1:lc)
      close(12)
      call copenfw(12, '.Revmax', icon)
      if(icon .ne. 0) then
         call cerrorMsg('.Revmax file cannot be opened', 0)
      endif
      write(msg, *) sngl( Zrevmax )
!        revmove blank at head ( gnuplot awk problem )
      call kseblk(msg, ' ', lc)      
      write(12,'(a)') msg(1:lc)
      close(12)
      end
!    *******************
      subroutine eplevel
      implicit none
#include "ZepDraw.h"
      integer i

 10   continue
       write(msg,*)  ' max level=',maxlevel,
     * ' current level range:  min=',clevelmin, ' max=',clevelmax
       call cerrorMsg(msg, 1)
       call cerrorMsg(' lvl:  numb. of comps', 1)
       do i = 0, min(maxlevel, maxmaxlevel)
          write(msg, '(i5, i8)') i, levelc(i)
          call cerrorMsg(msg, 1)
       enddo

       call cerrorMsg(
     * 'Enter min and max levels of comp. to be shown', 1)
       read(*,*) clevelmin, clevelmax
       if(clevelmin .lt. 0 .or. clevelmin .gt. maxlevel) goto 10
       if(clevelmax .lt. 0 ) goto 10
       if(clevelmax .lt. clevelmin) goto 10
       end
!    *********************
      subroutine epfixConfig
      implicit none
#include "ZepDraw.h"

      integer klena, icon, ios
      character*100 tempfile

      icon = 1

      tempfile = configFile
      do while (icon .ne. 0)
         write(0,'(a)')
     *   "Enter config file path: default is ",
     *   trim(configFile),
     *   "If above is OK, simply push Enter"
         read(*, '(a)',  iostat=ios) configFile
         if(configFile .eq. ' ') then
            configFile = tempfile
         endif
         
         if(ios .eq. 0) then
            call copenf(11, configFile, icon)
         else
            icon = 1   
         endif
         if(icon .ne. 0) then
            write(msg, *) configFile(1:klena(configFile))//
     *           ' cannot be opened; re-Enter the correct name '
            call cerrorMsg(msg, 1)
         else
            close(11)           ! it will  be reopened when it is read.
         endif
      enddo
      end

      subroutine epReadConfig
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
      integer i, j, w
       type(epPos)::  org
       type(ep3Vec)::  abc

!         read config file
      call eprcnf(configFile)

      call epqcnf(org, abc)
      write(0, *) 'Bounding box of the world'
      write(0, '(a,1p,3g15.7)') ' origing=', org
      write(0, '(a,1p,3g15.7)') ' box x y z length=',abc

      NumComp = Det%nct
      NumWorld = Det%nworld
      maxlevel = 0
      do i = 0, maxmaxlevel
         levelc(i) = 0
      enddo

      do i = 1, Det%nct
         j = min(maxmaxlevel, Det%cmp(i)%level)
         if(Det%cmp(i)%struc(1:3) .eq. 'box') then
            w = 1
         elseif(Det%cmp(i)%struc(1:6) .eq. 'sphere') then
            w = 20
         elseif(Det%cmp(i)%struc(1:3) .eq. 'cyl') then
            w = 10
         elseif(Det%cmp(i)%struc(1:4) .eq. 'pipe') then
            w = 20
         elseif(Det%cmp(i)%struc .eq. 'prism' .or.
     *          Det%cmp(i)%struc(1:6) .eq. 'prism_') then
            w = 2
         else
            w = 10
         endif
         levelc(j) = levelc(j) + w
         maxlevel = max(maxlevel, Det%cmp(i)%level)
      enddo
!      call mediaCount(Det.cmp, matdef, noOfMedia, compnToMn)
      call mediaCount(Det%cmp, Det%nct, matdef, noOfMedia)
      write(0,*) 'no of comps for each media'
      do i = 1, noOfMedia
         write(0,'(a,1x,i6)') matdef(i), NoOfCompsEachMedia(i)
      enddo
      end
      subroutine epMemoCnfgFile
!          memorize current config file name in Work/.configrc
      implicit none
#include "Zep3Vec.h"
#include "ZepDraw.h"
      integer icon
      call copenfw(21, '/tmp/$USER/.configrc', icon)
      if(icon  .eq. 0) then
         write( 21, '(a)' )  configFile
         close(21)
      endif
      end
      subroutine epOpenMediaFile
      implicit none
#include "ZepDraw.h"
      integer i, icon, klena
      character*8 mat
      integer,external:: kgetenv2
      character(24)::uname
!           open config media file:   config.pb etc 
      do i = 1, noOfMedia
         mat = matdef(i)
         icon = kgetenv2("USER", uname)
         fn(i) =
     *   '/tmp/'//uname(1:icon)//'/Work/'//'config.'//mat(1:klena(mat)) ! compose file name
         call copenfw(i + offset, fn(i), icon)   !  open it.
         if(icon .ne. 0) then
            open(i+ offset, file=fn(i), form='formatted')
         endif
         usedmedia(i) =.false.
      enddo
      end
!     ***********************
      subroutine epCloseDraw
      implicit none
#include "ZepDraw.h"
!
      integer i
      do i = 1, noOfMedia
         close(i+offset)
      enddo
      end

!     ***********************
      subroutine epexMedia(cmp, i, icon)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

       type(Component)::  cmp(*)
 
!          see if i-th component media is in the
!       selected list ?
! 
      integer i, icon
      integer k

      do k = 1, nselMedia
         if(seledMedia(k) .eq. cmp(i)%matter ) then
            icon = selmedia-1
            goto 10
         endif
      enddo
      icon =2 - selmedia
 10   continue
      end
      subroutine epexStruc(cmp, i, icon)
      implicit none
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"

       type(Component)::  cmp(*)

      integer i  ! input. i-th comp. is to be shown ?
      integer icon    ! output. 0. yes, 1 no.
      integer k

      do k = 1, nselstruc
         if(selstrucarray(k) .eq. cmp(i)%struc ) then
            icon = selstruc-1
            goto 10
         endif
      enddo
      icon =2 - selstruc
 10   continue
      end
!     ********************
      subroutine epbyStruc
      implicit none
#include "ZepDraw.h"
      
      integer klena, i, j

      call cerrorMsg('Available structures are as follows', 1)
      call cerrorMsg(struclist, 1)
      selstruclist = ' '

      
      selstruc = 0
      do while (selstruc .le. 0 .or. selstruc .gt. 2)
         call cerrorMsg(
     *   ' 1) Select structure to be shown',1)
         call cerrorMsg(
     *   ' 2) Select structure to be hidden',1)
         call cerrorMsg('Enter 1 or 2', 1)
         read(*,*) selstruc
      enddo
      write(msg, *)
     * "Enter structure name list, like: "//
     *  struclist(1:klena(struclist))//
     *  "  to show or to hide"
      call cerrorMsg(msg, 1)
      read(*, '(a)') selstruclist
      call kgetField(selstruclist, selstrucarray, maxsel, nselstruc)

      byStruc = .false.
      do i = 1, nselstruc
         do j = 1, nstruc
            if( strucarray(j) .eq. selstrucarray(i)) then
               byStruc =.true.
            endif
         enddo
      enddo

      if(.not. byStruc) then
         call cerrorMsg(
     *   'Your specification not exists in the config.', 1)
         call cerrorMsg('Your condition cancelled',1)
      endif
      end

