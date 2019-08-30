      subroutine cMCSconstForMedia(md, pm, aMCS)
!      for this media, sampling table is read and 
!      stored in aMCS.
!      It is assumed that pm = -1 comes first for a given md
!      and next call is for pm=1. So the file opened for
!      pm=-1 is kept open for pm=1. After it is read for
!     pm=1, the file will be closed.
      use modMCScontrol
!      use modXSecMedia
      use modMCS
      implicit none
#include "Zmedia.h"
#include "Zmanager.h"
#include "Zmanagerp.h"
      !  CosOrEpi is used from above
      type(epmedia),intent(in):: md  ! media 
      integer,intent(in):: pm   ! >0 for e+, < 0 for e-
      type(MCSconst),intent(out):: aMCS

      integer::icon, ie

      if( CosOrEpi == "epics") then
!     filename is declared in modTPXS
         if(MCSdir == " " ) then
            MCSdir ="$LIBLOFT/Data/MCS/"
         endif
!     filename="$EPICSTOP/Data/MCS/"//trim(md%name)
      else
         if(MCSdir == " " ) then
            MCSdir ="$LIBLOFT/Data/MCS/"
         endif
!        filename="$LIBLOFT/Data/MCS/"//trim(md%name)
      endif
      filename=trim(MCSdir)//trim(MCSparam)//"/"//trim(md%name)

      if( pm == -1 ) then
         call copenf(TempDev, filename, icon)
         if( icon /= 0 ) then
            write(0,*) trim(filename), ' cannot be opened for MCS'
            if( CosOrEpi == "epics") then
               write(0,*) 'If the path: '
               write(0,*)  trim(MCSdir)//trim(MCSparam)
               write(0,*)
     *         ' is correct, ',trim(md%name),' may be missing'  
               write(0,*) ' If so, you have to do next:'
               write(0,*) 
     *         '1: goto $LIBLOFT/Util/Elemag/MixedMCS'
               write(0,*) 
     *              '2: Edit paramdata file there: fix the 1st line'
               write(0,*)
     *         '   Normally only 1st term may be given; (0.05~1.9)'
               write(0,*) '3: then do'
               write(0,'(a,a,a,a)') '   ./ForManyMedia%sh ',
     *             ' $LIBLOFT/Data/MCS/',trim(MCSparam),
     *             " "//trim(md%name)
            endif
            stop
         endif
         read(TempDev, *)  ! skip 1 line  which is  C1forHardScat
         call creadMCSTab(TempDev, aMCS)
         read(TempDev, *)  
      elseif(pm == 1 ) then
         call creadMCSTab(TempDev, aMCS)
         close(TempDeV)
      else
         write(0,*)' pm =',pm, ' invalid for  cMCSconstForMedia'
         stop
      endif

      aMCS%loglambdah(:) =log( aMCS%lambdah(:) )

      ie = aMCS%minNon0mucEindex  ! for ie, 0 is stored
      aMCS%loglambdas1(ie:) =log( aMCS%lambdas1(ie:) )
      aMCS%loglambdas2(ie:) =log( aMCS%lambdas2(ie:) )
      aMCS%logmuc(ie:) =log( aMCS%muc(ie:) )

      end  subroutine cMCSconstForMedia
