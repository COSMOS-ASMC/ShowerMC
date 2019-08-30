! #include "ZsubstRec.h"
! #include "Zunionmap.h"
      subroutine epfixProc(id, media, path)
!     select the interaction process and the length in r.l or
!     kg/m2       
!          this is for EPICS       
      use modV1ry
      use modIntInf
      implicit none
#include  "Zglobalc.h"
#include  "Zmedia.h"
      
!#include  "Ztrack.h"
!     #include  "Ztrackv.h"
      integer,intent(in):: id   ! 1,2-> projetile is e,g
                                ! >2->  mu, pi,K, ...
      type(epmedia),intent(in):: media
      
      
!      real*8  den     ! input.  denstiy g/cm^3
      real(8),intent(out)::  path !  path length in r.l (EM=1,2)
!                      or in kg/m2  (EM> 2).
      
      real*8 len, minlen, den,  mTokgpm2
!     
      integer i
      
      minlen = Infty

      if( V1ry == 2 ) then
         if( id > 2 ) then
            call epForceV1ryInt !  make path for the target int.  0,
            !   so that next section select it.
         endif
      endif

      den = media%rho * media%rhoc    ! g/cm3    (g/1000)kg / (cm3/10^6) m^3 
      mTokgpm2 = 1.d3*den       ! m to kg/m2 conversion:
                                ! air  den-10^-3--> 1m=> 1kg/m2

      do i = 1, NumberOfInte
         if( IntInfArray(i)%decay ) then
            len = IntInfArray(i)%length
         else
            if( id <= 2) then
               len = IntInfArray(i)%thickness *  media%X0m / media%rhoc  !in m
            else
               len = IntInfArray(i)%thickness / mTokgpm2   ! in m
            endif
         endif
         if (len < minlen ) then
            ProcessNo = i
            minlen = len
         endif
      enddo

      path = minlen*1.d2
      i =   ProcessNo
      IntInfArray(i)%thickness  =  minlen * mTokgpm2 ! in kg/m2      
      IntInfArray(i)%length  =  minlen


      end
      
      subroutine epResetProcNoForV1ry
!          reset ProcessNo to be for the "hadint"  of virtual 1ry.
      use modV1ry
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
#include  "Ztrackv.h"
      integer:: i

      do i = 1, NumberOfInte

         if( IntInfArray(i)%process == hadint ) then 
            ProcessNo = i
            return  ! *******
         endif
      enddo
      write(0,*) ' for hadoron interaction=', hadint
      write(0,*) ' not in the interaction candidates in IntInfArray'
      write(0,*)
     *  ' Virtual 1ry treatment error: msg from  epResetProcNoForV1ry'
      stop
      end

      subroutine epsaveFirstCol
!        This is called only when the first nuclear interaction
!        or decay to save the collision info.
      implicit none
#include "Zepi2cos.h"      
#include "ZepMaxdef.h"


       type(ptcl)::  ptclSave(EPMAX_STACK)
      integer nptcls
      common /stackSavec/ ptclSave, nptcls
      
      integer i
      if( Nproduced >  EPMAX_STACK) then
         write(0,*) "Warning, first interacton generated: "
         write(0,*) "# of particls =",Nproduced, " > ",
     *               EPMAX_STACK, " so "
         write(0,*) " only first ", EPMAX_STACK, " ptcls are saved"
         write(0,*) 
     *     " for later retrieving (by calling epqFirstColPtcls)."
         write(0,*) " Simulation itself are not affected. "
         write(0,*)
     *     " If you need all the generated particle information "
         write(0,*) " please consider the use of epUI interface."
      endif
      nptcls = min(Nproduced, EPMAX_STACK)

      do i = 1,  nptcls
         ptclSave(i) = Pwork(i)
      enddo
      end
      subroutine epqFirstColPtcls(ptcls, n, m)
        implicit none
#include  "Zptcl.h"
#include "ZepMaxdef.h"
       type(ptcl)::  ptclSave(EPMAX_STACK)
      integer nptcls
      common /stackSavec/ ptclSave, nptcls
      integer n, m
       type(ptcl)::  ptcls(n)
      save
      if(n .lt. nptcls) then
         write(0,*)
     *    " n must be > ", nptcls, " in epqFirstColPtcls"
         stop  9999
      endif
      do m = 1, nptcls
         ptcls(m) = ptclSave(m)
      enddo
      m = nptcls
      end
