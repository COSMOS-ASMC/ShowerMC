          !This file contains stuff  for Light=21 or 22
          !
      subroutine epLightIOwriteIni
      implicit none
#include "ZsepManager.h"
          !     next is  needed for "Light"
#include "ZepTrackp.h"    
          ! used when Light is 21
      integer icon
      call epSeePMfile( OutPrimaryFile,  Out1ry)
      if( Out1ry > 0 ) then
         write(0,*) ' OutPrimaryFile=',trim(OutPrimaryFile)
         write(0,*) ' The file name must start with "-" for Light=21'
         write(0,*) ' Sorry. This may be redundant requirement but'
         write(0,*) ' is forced for future '
         stop
!         call copenfw2(OutPrimaryFileNo, OutPrimaryFile, 1, icon) ! ascii
      elseif( Out1ry < 0 ) then
         call copenfw2(OutPrimaryFileNo, OutPrimaryFile, 2, icon) ! binary
      else
         write(0,*) 'OutPrimaryFile name=',
     *    trim(OutPrimaryFile), 
     *   ' must have "-"(for binary) at its head'
         stop
      endif
      if(icon > 2 ) then
         write(0,*)  'Light=21, but ', trim(OutPrimaryFile)
         write(0,*)  ' could not be opened as file #=',
     *        OutPrimaryFileNo             
         stop
      endif
           ! open scratch file to store charged particle track
           ! info. for later Cerenkov generation.  This is used
           ! to store such tracks before the first interaction 
           ! occurs. (We want to put the first interaction
           ! info. in the OutPrimaryFile at its head)
      if(Out1ry > 0 ) then
!         open(IoScratch,  form="formatted",  status="scratch")
      else
         open(IoScratch,  form="unformatted",  status="scratch")
      endif
      

      if(Out1ry > 0 ) then
          ! write header info
 !        write(OutPrimaryFileNo,'(a)') 
 !     *     "#  mulsubTExyzdirwlmassCn"
 !         write(OutPrimaryFileNo,'(a)') 
 !     *        "#-------------------------------------" 
      else
         buf = "#  mulsubTExyzdirwlmassCn"
         write(OutPrimaryFileNo)  buf
         buf = "#-------------------"
         write(OutPrimaryFileNo) buf
      endif
      end
      
      subroutine epLightIOreadIni
      implicit none
#include  "ZsepManager.h"
#include  "ZepTrackp.h"

      integer icon

      if(Inp1ry == 1 ) then
         if( Light==22) then
            write(0,*) "PrimaryFile=",trim(PrimaryFile)
            write(0,*) "name must start with '-' for Light=22" 
            write(0,*) " sorry "
            stop
         endif
         call copenfw2(Ioprim, trim(PrimaryFile), 1, icon)
      elseif(Inp1ry == -1) then
         call copenfw2(Ioprim, trim(PrimaryFile), 2, icon)
      else
         write(0,*) ' file name and Inp1ry inconsistent'
         stop
      endif
      if(icon > 2 ) then
         write(0,*) ' +primry or -primary file=',trim(PrimaryFile)
         write(0,*) ' cannot be opened'
         stop
      endif
      if( Inp1ry > 0 ) then
         write(0,*) trim(PrimaryFile),
     *     ' is opened with formatted mode'
      else
         write(0,*) trim(PrimaryFile),
     *   ' is opened with binary mode'
      endif
      call epLightIOreadHD
      end

      subroutine epLightIOwriteiev
        ! init for each event
      implicit none
#include "ZsepManager.h"
      rewind IoScratch
              ! until first col. OutPrimEff is directed to IoScratch
      OutPrimEff = IoScratch
      end

      subroutine epwrite1stCol(cevent)
        ! Just after the first col. takes place 
      implicit none
#include "ZsepManager.h"
#include "ZepTrackv.h"
      integer,intent(in):: cevent

      integer::code, subcode, chg
      real(8):: erg, xin, yin,  zin
      real(8):: wx, wy, wz
      real(8):: wl, mass
      integer:: compno

      if(Out1ry > 0 ) then
!         write(OutPrimaryFileNo,*)  cevent
!         write(OutPrimaryFileNo,*) 
!     *     Incident.p.code, 
!     *     Incident.p.subcode, Incident.p.charge,
!     *     Incident.p.fm.p, Incident.p.mass,
!     *     Incident.pos,
!     *     Incident.t, Incident.w, Incident.wgt, Incident.pol,
!     *     Incident.wl, Incident.cn, Incident.user
!c/////////////
!         write(0,*) ' Incient=', 
!     *     Incident.p.code, 
!     *     Incident.p.subcode, Incident.p.charge,
!     *     Incident.p.fm.p, Incident.p.mass,
!     *     Incident.pos,
!     *     Incident.t, Incident.w, Incident.wgt, Incident.pol,
!     *     Incident.wl, Incident.cn, Incident.user
!c/////////
!         write(OutPrimaryFileNo,*)  FirstInt
!/////////////
!         write(0,*) ' FirstInt=',   FirstInt
!/////////
!         write(OutPrimaryFileNo,*)  FirstMedia.name, 
!     *    FirstMedia.colZ, FirstMedia.colA
!/////////////
!         write(0,*) ' FirstMedia=',   FirstMedia.name,
!     *    FirstMedia.colZ, FirstMedia.colA
!/////////
         write(OutPrimaryFileNo,*)  Proc1
!/////////////
!         write(0,*) ' Proc1=',    Proc1
!/////////
      else
         write(OutPrimaryFileNo)  cevent
         write(OutPrimaryFileNo)  Incident
         write(OutPrimaryFileNo)  FirstInt
         write(OutPrimaryFileNo)  FirstMedia%name, 
     *         FirstMedia%colZ, FirstMedia%colA
         write(OutPrimaryFileNo)  Proc1
      endif
      end
      subroutine epLightIOread1stCol(cevent, icon)
        ! when Light==22, get incident info.
      implicit none
#include "ZepManager.h"
#include "ZsepManager.h"
#include "ZepTrackv.h"
      integer,intent(out)::icon ! 0--> ok. -1-->end of all events (equiv.EOF)
      integer,intent(out):: cevent  ! event #
      if(Inp1ry > 0 ) then
         write(0,*) ' formatted primary is strange for Light=22'
         stop
!         read(Ioprim,*, end=100)  cevent
!         read(Ioprim,*) 
!     *     Incident.p.code, 
!     *     Incident.p.subcode, Incident.p.charge,
!     *     Incident.p.fm.p, Incident.p.mass,
!     *     Incident.pos,
!     *     Incident.t, Incident.w, Incident.wgt, Incident.pol,
!     *     Incident.wl, Incident.cn, Incident.user
!///////////
!         write(0,*)
!     *     Incident.p.code, 
!     *     Incident.p.subcode, Incident.p.charge,
!     *     Incident.p.fm.p, Incident.p.mass,
!     *     Incident.pos,
!     *     Incident.t, Incident.w, Incident.wgt, Incident.pol,
!     *     Incident.wl, Incident.cn, Incident.user
!////////////
!
         read(Ioprim,*)  FirstInt
         read(Ioprim,*)  FirstMedia%name,
     *         FirstMedia%colZ, FirstMedia%colA
!///////////
!         write(0,'(a,2f8.0)')  FirstMedia.name,
!     *         FirstMedia.colZ, FirstMedia.colA
!////////////

         read(Ioprim,*)  Proc1
!///////////
!         write(0,*)  Proc1
!////////////

      else
         read(Ioprim, end=100)  cevent
         read(Ioprim)  Incident
!///////////
!         write(0,*)
!     *     Incident.p.code, 
!     *     Incident.p.subcode, Incident.p.charge,
!     *     Incident.p.fm.p, Incident.p.mass,
!     *     Incident.pos,
!     *     Incident.t, Incident.w, Incident.wgt, Incident.pol,
!     *     Incident.wl, Incident.cn, Incident.user
!c////////////
         read(Ioprim)  FirstInt
         read(Ioprim)  FirstMedia%name,
     *         FirstMedia%colZ, FirstMedia%colA
!///////////
!         write(0,*) FirstMedia.name,
!     *         FirstMedia.colZ, FirstMedia.colA
!////////////
         read(Ioprim)  Proc1
      endif
      Nevrun = cevent
      if(cevent == 0 ) then
         icon = -1
      else
         icon = 0
      endif
      return
 100  continue
      icon = -1
      end
      

      subroutine epLightIOwritedE(io, nc)
      use modUI
            ! write energy deposit at every comp.
      implicit none
#include "ZsepManager.h"
      integer,intent(in):: io   !  file no.
      integer,intent(in):: nc   ! total comp. #
      integer i
      if( Out1ry > 0 ) then
         write(io,*) nc
      else
         write(io) nc
      endif
      do i=1,  nc
         if( ElossE(i) > 0. ) then
            if( Out1ry > 0  ) then
               write(io, *)  i, ElossE(i), ElossT(i)
            else
               write(io) i, ElossE(i), ElossT(i)
            endif
         endif
      enddo
      if( Out1ry >0 ) then
         write(io, *)  0,  0., 0.
      else
         write(io)  0,  0., 0.
      endif
      end

      subroutine epLightIOreaddE(totalcn, alloc)
      use modUI
            ! read energy deposit of each comp.
      implicit none
#include "ZsepManager.h"
      integer,intent(out) ::totalcn
            ! dirty workaround
      integer,intent(in) :: alloc  ! 0--> no allocation
                                       ! 1--> allcation
      integer::i
      real(4)::temp1, temp2  ! to avoid 0 index.
      if(Inp1ry > 0 ) then
         read(Ioprim,*,end=100) totalcn
      else
         read(Ioprim)  totalcn
      endif
      if( alloc > 0 ) then
         if( .not. allocated(ElossT)) then
            allocate(ElossT(totalcn) )
            allocate(ElossE(totalcn) )
         endif
      endif
      ElossT(:) = 0.
      ElossE(:) = 0.


      do while( .true. ) 
         if(Inp1ry > 0 ) then
            i = 0
!            read(Ioprim,*, end=100) i,  temp1, temp2
         else
            read(Ioprim, end=100) i,  temp1, temp2
         endif
         if(i == 0 ) exit
         ElossE(i) = temp1
         ElossT(i) = temp2
      enddo
 100  continue
      end


      subroutine epLightIOreadHD
          ! read  header part of the +primary or -primary
      implicitnone
#include  "ZsepManager.h"
#include  "ZepTrackp.h"

      integer icon, loc

      if(Inp1ry > 0 )  then
          read(Ioprim, '(a)')  buf
       else
          read(Ioprim)  buf
       endif

       if( index(buf, "mul") .gt. 0 ) then
          Inpmul = 1
       else
          Inpmul = 0
       endif
       if( index(buf, "xyz") .gt. 0  ) then
          Inpxyz = 1
       else
          Inpxyz = 0
       endif
       if( index(buf, "dir")  .gt. 0 ) then
          Inpdir =  1
       else
          Inpdir = 0
       endif
       if( index(buf, "sub") .gt. 0 ) then
          Inpsubcode = 1
       else
          Inpsubcode = 0
       endif
       if( index(buf, "KE") .gt.  0  .or.
     *     index(buf, "ke") .gt.  0   ) then
          Inperg = 1
       elseif( index(buf, "TE") .gt. 0 .or.
     *     index(buf, "te") .gt. 0  ) then
          Inperg = 2
       else
          Inperg = 0
       endif
       if(index(buf, "time") .gt. 0 ) then
          Inptime = 1
       else
          Inptime = 0
       endif

       if(index(buf, "Cn") .gt. 0 ) then
          Inpcn = 1
       else
          Inpcn = 0
       endif

       if(index(buf, "pol") .gt. 0 ) then
          Inppol = 1
       else
          Inppol = 0
       endif

       if(index(buf, "wgt") .gt. 0 ) then
          Inpwgt = 1
       else
          Inpwgt = 0
       endif

       if(index(buf, "wl") .gt. 0 ) then
          Inpwl = 1
       else
          Inpwl = 0
       endif

       if(index(buf, "mass") .gt. 0 ) then
          Inpmass = 1
       else
          Inpmass = 0
       endif

!       if(index(buf ,"Light") > 0 ) then
!          InpLight = 1
!       else
!          InpLight = 0
!       endif

       loc = index(buf, "user")
       if( loc .gt. 0) then
          Inpuser = 1
       else
          Inpuser = 0
       endif

       loc = index(buf, "disk")
       if( loc .gt. 0) then
          read(buf(loc+4:loc+5), *) Inpdisk
       else
          Inpdisk = 0
       endif

       write(0,*) '   Inpmul =', Inpmul 
       write(0,*) '   Inpxyz =', Inpxyz 
       write(0,*) '   Inpdir =', Inpdir
       write(0,*) 'Inpsubcode=', Inpsubcode
       write(0,*) '   Inperg =', Inperg
       write(0,*) '  Inptime =', Inptime 
       write(0,*) '    Inpcn =', Inpcn 
       write(0,*) '   Inppol =', Inppol 
       write(0,*) '   Inpwgt =', Inpwgt
       write(0,*) '    Inpwl =', Inpwl 
       write(0,*) '  Inpmass =', Inpmass 
!       write(0,*) ' InpLight =', InpLight
       write(0,*) '  Inpuser =', Inpuser 
       write(0,*) '  Inpdisk =', Inpdisk 
       call epskipcom(icon)     ! skip comment
       if(icon /= 0 ) then
          call cerrorMsg("Primary file has no #----- line", 0)
       endif
       end
