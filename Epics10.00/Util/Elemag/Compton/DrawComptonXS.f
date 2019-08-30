      implicit none
!
!         total prob. of compton effect 
!
#include "Zglobalc.h"
#include "ZbasicCnst.h"      
#include "ZepManager.h"
#include "Zmedia.h"
#include "Zmass.h"
#include "ZepTrackp.h"
! #include "ZBPgene.h"

      character*130  file
      character(len=8)::name
      integer io, result

      real*8  Eg, E1, E2, step, xprob(5), txray(5)
      real*8  tprob, p
       type(epmedia)::  media
      integer XCOM, icon

      io = 10

      E1=1.e-6
 1    E2=100.
      step = 0.05
      write(0,*) 'Enter  XCOM(0/1) and media name'
      write(0,*) '   If XCOM=1, use XCOM xsec'
      write(0,*) '   If XCOM=0, not use XCOM xsec'
      read( *, * )  XCOM, name
      
      file = "$EPICSTOP/Data/Media/"//trim(name)
      write(0,*) 'file=',trim(file)

      call copenf(io, file, result)
      if( result .eq. 0 ) then
         call epReadMTbl(io, media)
         call epGetEffZA(media)
         MediaDir(1) ="$EPICSTOP/Data/Media"
         if(XCOM .eq. 1) then
            Excom1 = 100.
            call epReadXXsec(media, icon)
            if(icon .ne. 0) then
               write(0,*) ' Xcom read err'
            endif
         elseif(XCOM .eq. 0) then
            Excom1 = E1
         else
            write(0,*) ' XCOM is invalid =',XCOM
            stop 1234
         endif

         write(0,*) "#  'prob. of Compton'"
         write(0,*)
     *    "Eg(MeV),    /rl    /(g/cm2)    /cm     mb "
         Eg = E1
         do while ( Eg .lt. E2)   
            if(XCOM .eq. 0) then
               call epcompp(media, Eg, tprob, p) ! tprob is prob/r.l
            elseif( media%xcom%size .gt. 0 .and.
     *              Eg .le. Excom1 ) then
               call epXrayp(media, Eg, 1, 3,  xprob, txray)
               tprob= xprob(2)  ! /r.l
            else
               write(0,*) "media%xcom%size=",media%xcom%size
               write(0,*) " err in XS%out"
               stop
            endif
               write(*, '(1p,5g12.4)') 
     *       Eg*1000.,  
     *       tprob, tprob/media%X0g,
     *       tprob/media%X0, tprob/media%mbtoPX0
             Eg = Eg * 10.0d0**step
         enddo
      else
         write(0,*) ' media file error: file= ',trim(file)
      endif
      end
