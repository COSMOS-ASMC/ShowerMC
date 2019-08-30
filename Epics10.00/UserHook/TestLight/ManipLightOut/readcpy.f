#include "ZepicsBD.h"
      module modepreadcpy
      integer:: narg, leng

      integer::  code, subcode, chg, compno
      real(8):: erg,  xin, yin, zin,  wx, wy, wz,
     *           wl,  mass
      character(len=120)::evlistfile
      integer::totalcn, eventn
      logical:: copy, evlist
      integer,save::listfileno=9, fno
      integer::lineC     ! Cerenkov outuput lines
      integer::lineS         ! Scinti //
      integer::current

      end module modepreadcpy

      program main
      use modepreadcpy
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
#include  "ZsepManager.h"
      integer::icon
      integer:: i      
      call procCmdLine

      call epLightIOreadHD
      if(copy) then
         call epwriteHD
      endif
      current = 0
      do while (.true. )
         call epLightIOread1stCol(eventn, icon)
         if(icon /= 0 ) exit
         if( evlist ) then
            do while(  current < eventn) 
               read(listfileno, *, end=100) current
            enddo
         else
            current = eventn
         endif
         if(current == eventn )  then
            call epPrint1stCol(eventn)
            if(copy) then
               call epwrite1stCol(eventn)
            endif
         endif
         call procCerenScin(icon)

         call epLightIOreaddE(totalcn, 1) ! 1-->alloc ElossE,ElossT, if needed
         if( current == eventn ) then
            if(copy) then
               call epLightIOwritedE(OutPrimaryFileNo, totalcn)
            endif
            call epPrintEloss
         endif
      enddo
 100  continue
      end
      subroutine outmessg
      implicit none
      character(len=70)::msg(12)=(/
     *"Usage: ./read.. PrimaryFile [OutPrimaryFile eventlist detail]",
     *'PrimaryFile: output file from Light=21',
     *'OutPrimaryFile: File with "-" at its top to store copy of ',
     *'   selected events; If "dummy" or not given, no copy is made',
     *'eventlist:', 
     *'   file name which has a list of event # (one in one line)',
     *'   to be processed (selection/printing); if dummy or not given',
     *'   all events are target ',
     *'detail: if not given, info for charged particle track and ',
     *'   energy deposit  in a cell is printed only limited number',
     *'   else if a number is given, max of such number of particles',
     *'   is printed. If 0 is given, all particles will be printed '
     * /)
      write(0,*) msg
      end
      subroutine epPrint1stCol(eventn)
      implicit none
#include  "ZepTrackv.h"
c#include  "ZepManager.h"
c#include  "ZsepManager.h"

      integer,intent(in)::eventn
      write(*,*) eventn
      write(*,*) Incident.p.code, Incident.p.subcode, Incident.p.charge
      write(*,*) Incident.p.fm.p(4)
      write(*,*) Incident.pos
      write(*,*) Incident.w
      write(*,*) FirstInt
      write(*,*) FirstMedia.name,  FirstMedia.colZ, FirstMedia.colA
      write(*,*) Proc1

      end
      subroutine procCmdLine
      use modepreadcpy
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
#include  "ZsepManager.h"
      integer::icon
      character(len=20)::detail
      
      narg=iargc()
      if( narg < 1  .or. narg > 4 ) then
         call outmessg
         stop
      endif
      call getarg(1, PrimaryFile, leng)
      write(0,*) leng, trim(PrimaryFile)
      call epSeePMfile(PrimaryFile, Inp1ry)
      call copenfw2(Ioprim, trim(PrimaryFile), 2, icon)
      if( narg > 1 ) then
         call getarg(2, OutPrimaryFile, leng)
         write(0,*) leng, trim(OutPrimaryFile)
         copy =  trim(OutPrimaryFile) /= 'dummy'
         if(copy) then
            call epSeePMfile( OutPrimaryFile,  Out1ry)
            call copenfw2
     *          (OutPrimaryFileNo, trim(OUtPrimaryFile), 2, icon)
         endif
      else
         copy = .false.
      endif
      
      if(narg > 2 ) then
         call getarg(3, evlistfile, leng)
         if( trim(evlistfile) == 'dummy' ) then
            evlist = .false.
         else
            evlist = .true.
            call copenfw2(listfileno, evlistfile, 1, icon)
         endif
      else
         evlist = .false.
      endif
      if(narg == 4) then 
         call getarg(4, detail, leng)
         read(detail, *) fno
         if(fno == 0) then
            fno = 100000000
         endif
      else
         fno= 10
      endif
      end
      
      subroutine procCerenScin(icon)
      use modepreadcpy
      implicit none
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
#include  "ZepManager.h"
#include  "ZsepManager.h"

      integer::icon

      lineC = 0                ! Cerenkov outuput lines
      lineS = 0                ! Scinti //
      do while(.true.) 
         read(Ioprim,  end = 100)  
     *        code,
     *        subcode, chg, erg,  xin, yin, zin,
     *        wx, wy, wz,
     *        wl, mass,  compno
         if(current == eventn) then
            if(copy) then
               write(OutPrimaryFileNo) 
     *              code,
     *              subcode, chg, erg,  xin, yin, zin,
     *              wx, wy, wz,
     *              wl, mass, compno
            endif
            if(code == -3 ) then
               lineC = lineC + 1
               if(lineC <= fno ) then
                  write(*, 
     *                '(i6,2i3, 1p, 4g13.5, 3g15.8, 2g12.4, i6)')
     *                 code,
     *                subcode, chg, erg,  xin, yin, zin,
     *                wx, wy, wz,
     *                wl, mass, compno
                  if(lineC == fno)  then
                     write(*,*)
                  endif
               endif
            elseif( code == -2 ) then
               lineS = lineS + 1
               if(lineS <= fno ) then
                  write(*, 
     *                '(i6,2i3, 1p, 4g13.5, 3g15.8, 2g12.4, i6)')
     *                 code,
     *                subcode, chg, erg,  xin, yin, zin,
     *                wx, wy, wz,
     *                wl, mass, compno
                  if(lineS == fno)  then
                     write(*,*)
                  endif
               endif
            endif
         endif
         if(code == -1000)  exit
      enddo            
      icon = 0
      return
 100  continue
      icon  = -1
      end
      subroutine epPrintEloss
      use modepreadcpy
      implicit none
      real(8):: dEeff, dEt
      integer:: i
      integer:: counter
      
      counter = 0
      do i = 1, totalcn
         call epqEloss(i, dEeff, dEt)
         if(i == 0 ) exit
         if(dEeff > 0 ) then
            counter = counter + 1
            if(counter <= fno ) then
               write(*,*) i, dEeff, dEt
            endif
         endif
      enddo
      end

      subroutine epwriteHD
      implicit none
#include  "ZepTrackp.h"
c#include  "ZepManager.h"
#include  "ZsepManager.h"
      buf = "#  mulsubTExyzdirwlmassCn"
      write(OutPrimaryFileNo)  buf
      buf = "#-------------------"
      write(OutPrimaryFileNo) buf
      end

                                    
