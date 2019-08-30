      subroutine epOutCnfWithQuench(l1, l2, out)
      use epModify
       implicit none
#include  "ZepTrackp.h"
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
! #include  "Zptcl.h"
#include  "Zmass.h"
!#include  "Zmedia.h"
!         
       type(epmedia)::  mediax 
       integer l1, l2 ! input. level min, max
       integer out  ! input. output logical device no.
       integer i, j, m, k
       integer  klena
       character*200 f1, f2, f3, f4
       character*700 f5
       character*1024 temp
       integer lvlskip, jmax, j1 
       integer loc(maxiattr)
       character*2 bq/'\\'/  ! solaris cannout use '\'
       character*20 matterAndrhoc
       real rhoc
       integer lc, modf
       logical needpr
       real*8 dedx, c1, c2, cc
       character*1 id
       type(ptcl)::  aPtcl   ! dummyn
       dedx = 0.  ! dummy
       
       lvlskip = 0
       do  i = 1,  Det%nct
          modf = Det%cmp(i)%modifier
          MediaNo = Det%Cn2media(i)

          if(modf > 0 .or. Media(MediaNo)%Birks /=' ') then
             mediax = Media(MediaNo)
             call epQuenchCoeff(modf, mediax, aPtcl, dedx, c1,c2,cc,id)
             if(id == 'n') cycle
          else
             cycle
          endif
          if(wrtcom .ne. 0) then
             call epwrtcnfcom(out, i)
          endif
          if( ( Det%cmp(i)%level .ge. l1) .and. 
     *        ( Det%cmp(i)%level .le. l2) ) then
             if(lvlskip .ne. 0) then
                write(out, *)
     *          lvlskip, ' comps are skipped: those ' //
     *          ' not within the levels, (', l1, '~',l2,')'
                lvlskip = 0
             endif

             call epqmatrhoc(i, matterAndrhoc, lc, rhoc)
             write(temp,
     *         '(i6, a1, i3, a1, a12, a1, a, i2, i6)')  
     *         i, ' ', Det%cmp(i)%level, ' ', trim(Det%cmp(i)%struc),
     *         ' ',  matterAndrhoc(1:lc), Det%cmp(i)%CountIO, 
     *         Det%cmp(i)%CountDE
!             compress if 2.030000E-5 etc to 2.03E-5
             call epsup00E(temp)

             f1 = ' '
             call epcompBlank(temp, f1)

             temp = ' '
             if(  Det%cmp(i)%MaxPathL .le. 0.d0 .and.
     *            Det%cmp(i)%Modifier .eq. 0
     *              .and. Det%cmp(i)%chno .eq. NotGiven ) then
                write(temp, '(a6,3g27.10)') ' 0 0 /', 
     *               Det%cmp(i)%orgx,  Det%cmp(i)%orgy,
     *               Det%cmp(i)%orgz
             elseif( Det%cmp(i)%MaxPathL .gt. 0.d0 .and.
     *            Det%cmp(i)%Modifier .eq. 0
     *              .and. Det%cmp(i)%chno .eq. NotGiven ) then
                write(temp, '(1p,g15.5, a4,3g27.10)')
     *               Det%cmp(i)%MaxPathL,  ' 0 /', 
     *               Det%cmp(i)%orgx,  Det%cmp(i)%orgy,
     *               Det%cmp(i)%orgz
             elseif( Det%cmp(i)%MaxPathL .le. 0.d0 .and.
     *            Det%cmp(i)%Modifier > 0
     *              .and. Det%cmp(i)%chno .eq. NotGiven ) then
                write(temp, '(a3,i6, a3,1p, 3g27.10)')
     *                 ' 0 ',   Det%cmp(i)%Modifier," / ",
     *               Det%cmp(i)%orgx,  Det%cmp(i)%orgy,
     *               Det%cmp(i)%orgz

             elseif( Det%cmp(i)%chno .eq. NotGiven) then
                write(temp, '(1p,g15.5, i6, a2,3g27.10)')
     *               Det%cmp(i)%MaxPathL, Det%cmp(i)%Modifier,
     *               ' /', 
     *               Det%cmp(i)%orgx,  Det%cmp(i)%orgy,
     *               Det%cmp(i)%orgz
             elseif( Det%cmp(i)%chno .ne. NotGiven) then
                if( Det%cmp(i)%MaxPathL .le. 0.) then
                   write(temp, '(a3,i6, i9, a2, 3g27.10)')
     *                  ' 0 ',    Det%cmp(i)%Modifier,
     *                  Det%cmp(i)%chno, ' /',    
     *                  Det%cmp(i)%orgx,  Det%cmp(i)%orgy,
     *                  Det%cmp(i)%orgz
                else
                   write(temp, '(g15.5,i6, i8, a2, 3g27.10)')
     *                  Det%cmp(i)%MaxPathL, Det%cmp(i)%Modifier,
     *                  Det%cmp(i)%chno,  ' /',
     *                  Det%cmp(i)%orgx,  Det%cmp(i)%orgy,
     *                  Det%cmp(i)%orgz
                endif
             endif
             call epsup00E(temp)
             f2 = ' '
             call epcompZero(temp, f2)

             if( Det%cmp(i)%strucNo .gt. 0) then
                k = Det%cmp(i)%strucNo
                temp = ' '
                call epatloc(Det%cmp(i), loc)
                write(temp, '(50G14.5)')
     *               (Volat( Det%cmp(i)%vol+loc(j) ), j=1, 
     *                Det%cmp(i)%Nattributes)
!     *    or         nattrib(k))
             else
!                should be sub Detector
                write(temp, *) ' '
             endif
             call epsup00E(temp)
             f3 = ' '
             call epcompZero(temp, f3)

             if( .not. Det%cmp(i)%rotation .or.
     *            ( Det%cmp(i)%direc(1)  .eq. 1.d0 .and.
     *            Det%cmp(i)%direc(5)  .eq. 1.d0 .and.
     *            Det%cmp(i)%direc(9)  .eq. 1.d0) ) then
                write(temp, *) ' '
             else
                temp=' '
                write(temp,'(9G27.15)')
     *               ( Det%cmp(i)%direc(j),j=1,9)
             endif
             call epsup00E(temp)
             f4 = ' '
             call epcompZero(temp, f4)

             j = Det%cmp(i)%NMatreska
             if(j .gt. 0)  then
                jmax = min(j, 140)
                temp=' '
                write(temp, '(a2,1x,140i7)')
!((((((((((((((((
!     *        ' /', (Det.cmp(i).ContainsR(m),m=1,j)
     *           ' /', (CnArea( Det%cmp(i)%ContainsR+m ),m=1,jmax)
!))))))))))))))
                if(jmax .lt. j) then
!                   put \ at the end
                   temp(klena(temp)+2:klena(temp)+2)=bq(1:1)  ! continues
                endif
                f5 = ' '
                call epcompBlank(temp, f5)
             else
                jmax = 0
                f5 =' '
             endif
             write(out, *) f1(1:max(klena(f1),1)), 
     *            f2(1:max(1,klena(f2))), 
     *            f3(1:max(1,klena(f3))), 
     *            f4(1:max(1,klena(f4))),
     *            f5(1:max(1,klena(f5)))

             if(jmax .lt. j) then
!                 put remaining matreshka
                j1 = jmax + 1
                do while(j1 .le. j)
                   jmax = min(j1+140, j)
                   temp=' '
                   write(temp, *) 
     *             (CnArea( Det%cmp(i)%ContainsR+m ),m=j1,jmax)
                   if(jmax .lt. j) then
!                     put \ at the end
                      temp(klena(temp)+2:klena(temp)+2)=
     *                     bq(1:1)   ! continues
                   endif
                   write(out,'(a)') temp(1:klena(temp))
                   j1 = jmax+1
                enddo
             endif
!             call epSetEmin(i)
!             write(out,'("# min KE(keV):", 5f10.1)')
!     *       1.d6*Det.cmp(i).EminG, 1.d6*(Det.cmp(i).EminE-masele), 
!     *       1.d6*Det.cmp(i).RecoilE, 1.d6*KEmin, 1.d6*EminH
             if(id =="L") then
                write(out,'("# quench info: ",a, 1p, 3g12.3)')
     *             id, c1, c2, cc
             elseif(id == "T") then
                write(out,'("# quench info: ",a, 1p, 2g12.3)')
     *             id, c1, c2
             elseif(id == "B") then
                write(out,'("# quench info: ",a, 1p, 1g12.3)')
     *             id, c1
             endif
          else
!             skip
             lvlskip = lvlskip + 1
          endif
       enddo

       end
      subroutine epsuplast0(text)
      implicit none

!         supress last 0's such as in 0.0900 or 1.1000
!         and get  0.09 or 1.1
!      
      character*(*) text   ! in/out
      integer klena
      integer i
      i = index(text, ".")
      if( i .gt. 0 )   then
         i = klena(text)
         do while( text(i:i) .eq. "0") 
            text(i:i)=" "
            i = i - 1
         enddo
      endif
      end

!c        test epsup00E
!      character*100 text
!      text =
!     *   '0.50000000000E-02  0.5500000e1 1.00300e-2 1.33000000E-02'
!      call epsup00E(text)
!      write(*,*) text
!      end
      subroutine epsup00E(text)
      implicit none
!        when a config file output is made by epOutCnf,
!        some data may be represented as 
!        5.0000000000000E-1 which  can be shortned without
!        loss of accuracy as
!        5.E-1.
!        This routine shortns the string in text, if it has
!        0E  or 0e to (E or e). This process is repeated
!        until there is no 0E or 0e.
!
      character*(*) text  ! in/out.
      integer klena
      integer i, tl

      i = 1
      tl = klena(text)
      do while (i .ne. 0  .and. tl .gt. 7)
         i = index(text(1:tl), '0000000E')
         if(i .gt. 0) then
            text(i:tl) = text(i+7:tl)
         else
            i = index(text(1:tl), '0000000e')
            if(i .gt. 0) then
               text(i:tl) = text(i+7:tl)
            endif
         endif
         tl = klena(text)
      enddo

      i = 1
      do while (i .ne. 0 .and. tl .gt. 3) 
         i = index(text(1:tl), '00E')
         if(i .gt. 0) then
            text(i:tl) = text(i+2:tl)
         else
            i = index(text(1:tl), '00e')
            if(i .gt. 0) then
               text(i:tl) = text(i+2: tl)
            endif
         endif
         tl = klena(text)
      enddo


      i = 1
      do while (i .ne. 0 .and. tl .gt. 2) 
         i = index(text(1:tl), '0E')
         if(i .gt. 0) then
            text(i:tl) = text(i+1:tl)
         else
            i = index(text(1:tl), '0e')
            if(i .gt. 0) then
               text(i:tl) = text(i+1: tl)
            endif
         endif
         tl = klena(text)
      enddo
         
      end

      subroutine epcompBlank(x, y)
!       compress two or  more blanks to one blank.
!       y must be put blank beforehand.
      implicit none
      character*(*) x, y
      integer klena, n, i, nb, j
!         if y is x, the tail part is not
!         well treated.
!
      n = klena(x)
      nb = 0
      j = 0
      do i = 1, n
         if(x(i:i) .eq. ' ') then
            if(nb .eq. 0) then
               nb = nb + 1
               j = j + 1
               y(j:j) = x(i:i)
            else
               goto 10
            endif
         else
            j = j + 1
            y(j:j) = x(i:i)
            nb = 0
         endif
 10      continue
      enddo
      end
!     ***************
      subroutine epwrtcnfcom(io, i)
      implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
      integer io ! input.  to which comment is directed. (log. dev)
      integer i  ! input.  component number to be printed nextepr

      integer j
      
      do j = comloc, comcounter
         if(comflag(j) .gt. i) goto 10
         if(comflag(j) .eq. i)  then
            write(io,'(a)') comarea(j)
         endif
      enddo
 10   continue
      comloc = j
      end
!     ************************
      subroutine epatloc(comp, loc)
       implicit none
#include  "Zep3Vec.h"
#include  "Zcnfig.h"
!         

       type(Component)::  comp
       integer loc(*)

       integer i
       integer:: uscl

       character*8 epparaphrase, tempph
       character*80 msg

       if(comp%struc(1:3) .eq. 'box' .or.
     *    comp%struc(1:3) .eq. 'cyl' .or.          
     *    comp%struc(1:4) .eq. 'pipe' .or.          
     *    comp%struc .eq. 'prism' .or.          
     *    comp%struc(1:6) .eq. 'prism_' .or.          
     *    comp%struc(1:6) .eq. 'sphere' ) then
          do i = 1, comp%Nattributes
             loc(i) = i
          enddo
       else
          call epseeUnderScore(comp%struc, uscl)
          tempph = epparaphrase(comp%struc(1:uscl))
          if(tempph(1:4) .eq. 'new-') then
             call epatlocNew(comp, loc)
          else
             write(msg,*)
     *            'struc=', comp%struc,' not supported.'
             call cerrorMsg(msg,0)
          endif
       endif
      end
