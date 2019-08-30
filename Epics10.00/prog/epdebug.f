      subroutine epStrange
      implicit none
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"
#include "Zcode.h"

      integer::k
      integer,save::nstrange = 0
      nstrange = nstrange + 1                  

      k = cTrack%p%code
      write(0,*) '|Z|>1 but not heavy: code,sub,chg='
      write(0,*) k, cTrack%p%subcode, cTrack%p%charge
      write(0,*) 'K%E=',cTrack%p%fm%p(4)-cTrack%p%mass
      write(0,*) ' mass=',cTrack%p%mass
      cTrack%p%fm%p(4)=cTrack%p%mass  ! make 0 energy
      if(nstrange > 100 ) then
         stop
      endif
      write(0,*) nstrange, ' times neglected'
      end
      subroutine debugm(msg)
      implicit none
#include "ZepTrackv.h"
      character*(*) msg
      logical strange
      character*16 trigger

      strange = .false.
      if( Media(MediaNo)%noOfElem .le. 0 .or. 
     *    Media(MediaNo)%noOfElem .gt. 5 ) then
         strange = .true.
         trigger='noOfElem'
      elseif(Media(MediaNo)%A .lt. 1. .or.
     *   Media(MediaNo)%A .gt. 1000.) then
         strange = .true.
         trigger='A'
      elseif(Media(MediaNo)%Z .lt. 1. .or.
     *   Media(MediaNo)%Z .gt. 1000.) then
         strange = .true.
         trigger='Z'
      elseif(Media(MediaNo)%Aeff .lt.  1.) then
         strange = .true.
         trigger='Aeff'
      endif
!
      if(strange) then
         if( Media(MediaNo)%rhoc .ne.  1.) then
            write(*,*) 'No rhoc correction included in the next'
         endif
         write(*,*) trigger, ' ***********', msg
         write(*,*) ' code=',cTrack%p%code,' proc=',Move%proc
         write(*,*) '      MediaNo=', MediaNo,
     *      ' noOfElem=',Media(MediaNo)%noOfElem
         write(*,*)
     *         ' charge=',cTrack%p%charge, ' subcode=',cTrack%p%subcode,
     *         ' E=', cTrack%p%fm%p(4), ' Zeff=', Media(MediaNo)%Zeff
         write(*, *) ' noOfElem=', Media(MediaNo)%noOfElem,
     *         ' n, A, Z', Media(MediaNo)%n,
     *          Media(MediaNo)%A, Media(MediaNo)%Z
         write(*,*)' mbtoPgrm=',Media(MediaNo)%mbtoPgrm,
     *         ' mbtoPcm=',Media(MediaNo)%mbtoPcm,
     *         ' mbtoPX0=',Media(MediaNo)%mbtoPX0

         write(*,*)' mbtoPgrm2=',Media(MediaNo)%mbtoPgrm2,
     *         ' mbtoPcm2=',Media(MediaNo)%mbtoPcm2,
     *          ' mbtoPX02=',Media(MediaNo)%mbtoPX02

         write(*,*)' Z2byAeff, Z5byA=', Media(MediaNo)%Z2byAeff,
     *         Media(MediaNo)%Z5byAeff 


         write(*,*)  ' A=', Media(MediaNo)%A,
     *           ' Z=', Media(MediaNo)%Z,
     *           ' Aeff=', Media(MediaNo)%Aeff
      endif
      end
      subroutine  debugpos(msg, compn, pos, icon)
      implicit none
#include "ZepTrackv.h"
       type(epPos)::  pos
      character*(*) msg
      integer compn, icon
         write(*,*) '--------------- ', msg, ' Cn=',Cn
         write(*,*) ' compn=',compn, ' pos%x,y,z=',
     *      pos%x, pos%y,  pos%z
         write(*,*) ' ctrackpos=',cTrack%pos%x, 
     *      cTrack%pos%y, cTrack%pos%z, ' wx,y,z=',
     *      cTrack%w%x, cTrack%w%y, cTrack%w%z
         write(*,*) ' M-w%x =',Move%Track%w%x, 
     *        Move%Track%w%y, Move%Track%w%z
         write(*,*) ' Move%track%pos=',Move%Track%pos%x,
     *     Move%Track%pos%y,  Move%Track%pos%z, ' dl=',
     *     Move%dl, Move%dE, Move%Track%p%mass,
     *    ' trucn=',Move%Trunc
         write(*,*)
     *     ' code=', cTrack%p%code,  ' chg=',cTrack%p%charge,
     *     ' E=',  cTrack%p%fm%p(4)
         write(*,*) ' cross=',Move%Cross, ' Trunc=',Move%Trunc
         write(*,*) ' icon=',icon
      end
      subroutine  epfordebug(msg)
      implicit none
#include "ZepTrackv.h"
       type(epPos):: postemp
       type(epDirec)::  dirtemp
      character*(*) msg  ! callers name is usually better
      character(len=20) struc

         write(0,*) '--------------- ', msg
         write(0,*) ' ctrackpos=',cTrack%pos%x, 
     *      cTrack%pos%y, cTrack%pos%z, ' wx,y,z=',
     *      cTrack%w%x, cTrack%w%y, cTrack%w%z
         write(0,*) ' M-w%x =',Move%Track%w%x, 
     *        Move%Track%w%y, Move%Track%w%z
         write(0,*) ' Move%track%pos=',Move%Track%pos%x,
     *     Move%Track%pos%y,  Move%Track%pos%z, ' dl=',
     *     Move%dl,' dE and mass=', Move%dE, Move%Track%p%mass,
     *    ' trucn=',Move%Trunc
         write(0,*)
     *     ' code=', cTrack%p%code,  ' chg=',cTrack%p%charge,
     *     ' E=',  cTrack%p%fm%p(4), ' wl =',cTrack%wl
         write(0,*) ' cross=',Move%Cross, ' Trunc=',Move%Trunc
         write(0,*) ' cTrack%cn=',cTrack%cn, ' Move%cn=',
     *    Move%track%cn, 'Cn=',Cn
         write(0,*) ' media #=', MediaNo, ' media=', Media(MediaNo)%name
         call epqstruc(Cn, struc)
         write(0,*) ' struc =', struc
         write(0,*) ' Move%proc=', Move%proc
!             world coord.

         call epl2w(Cn, cTrack%pos, postemp)
         call epl2wd(Cn, cTrack%w,  dirtemp)
         write(0,*) ' assuming Cn=',Cn, ' world pos=',
     *   postemp%x, postemp%y, postemp%z
         write(0,*) ' dir=', dirtemp%x, dirtemp%y, dirtemp%z
         write(0,*)
      end
      subroutine epcurrent
#include "ZepTrackp.h"
#include "ZepTrackv.h"
            write(*,*) ' now ground'
            write(*,*) ' target A,Z=',Media(MediaNo)%colA,
     *         Media(MediaNo)%colZ
            write(*,*) ' code=',cTrack%p%code, ' chg=',
     *                 cTrack%p%charge, ' sub=', cTrack%p%subcode,
     *               ' Ek=',cTrack%p%fm%p(4)-cTrack%p%mass
            end
