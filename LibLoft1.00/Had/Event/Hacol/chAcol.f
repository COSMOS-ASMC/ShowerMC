!         hadoron-A collision.  generic.
      subroutine chAcol(pj, ia, iz, xs, a, ntp)
      implicit none
#include  "Zcondc.h"
#include  "Zcode.h"
#include  "Zptcl.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zspecial.h"
!
      type(ptcl):: pj   ! input  projectile ptcl
      integer ia    ! input. mass no. of target
      integer iz    ! input. charge no. of target
      real(8),intent(in):: xs ! x-section for nucleus ia (mb)
      type(ptcl):: a(*)  !  output. produced ptcls
      integer ntp   ! number of produced ptcls

      integer nc
      integer retry
      data nc/0/
      save nc
      retry = 1
      do while( retry .ne. 0) 
         call chAcolx(pj, ia, iz, xs, a, ntp)
         if(modifyx .ne. 0) then
            call ccheckx(pj, ia, iz, a, ntp, retry)
            nc = nc + 1
            if(nc .lt. 20) then
               write(0,*)
     *         '************* X modification is requested !!!!!!!!!'
            endif
         else
            retry = 0
         endif
      enddo
      end
      subroutine ccheckx(pj, ia, iz, a, ntp, retry)
      implicit none
#include  "Zcondc.h"
#include  "Zcode.h"
#include  "Zptcl.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zspecial.h"
!
      type(ptcl):: pj   ! input  projectile ptcl
      integer ia    ! input. mass no. of target
      integer iz    ! input. charge no. of target
      type(ptcl):: a(*)  !  output. produced ptcls
      integer ntp   ! number of produced ptcls
      integer retry  !  output. 0 if event is acceptable

      integer i
      real*8 u, ePimax, xPimax, eNmax, xNmax, xxx

      if(pj%code .eq. knuc) then
!         find max nucleon
         ePimax = 0.
         eNmax = 0.
         do i = 1, ntp
            if( a(i)%code .eq. knuc  )then
               eNmax = max(eNmax, a(i)%fm%p(4) )
            elseif( a(i)%code .eq. kpion  )then
               ePimax = max(ePimax, a(i)%fm%p(4) )
            endif
         enddo
         xPimax = ePimax/pj%fm%p(4)
         xNmax = eNmax/pj%fm%p(4)
         call rndc(u)
         if(xNmax .gt. 0.97) then
            retry = 0
         elseif(modifyx .eq. 1) then
!            drop  high Pi x; pw=3
            if(u .lt. (1.0 - xPimax)**modifyxpw1) then
               retry = 0
            else
               retry = 1
            endif
         elseif(modifyx .eq. 2) then
!              drop high Nx case; pw=3
            if( u .lt. ((2.-xNmax)/2.)**modifyxpw1 ) then
               retry = 0
            else
               retry = 1
            endif
         elseif(modifyx .eq. 3) then
!               drop large xPimax*xNmax ; pw =3
            xxx = xPimax* xNmax
            if(xxx .lt. 0.005) then
               retry = 0
            elseif( xxx .gt. 0.08 ) then
               retry = 1
            elseif(u .lt. ((0.08-xxx)/(0.08-0.005))**modifyxpw1 ) then
               retry = 0
            else
               retry = 1
            endif
         elseif(modifyx .eq. 4) then
            xxx = xPimax/xNmax
            if(xxx .gt. 2.5) then
               retry = 1
            elseif(xxx .lt. 0.2) then
               retry = 0
            elseif(u .lt. ((2.5-xxx)/(2.5-0.2))**modifyxpw1) then
               retry = 0
            else
               retry = 1
            endif
         elseif(modifyx .eq. 5) then
            xxx=xPimax/xNmax
            if(xxx .gt. 2.)  then
               retry = 0
            elseif( u .gt.  (1.-xxx/2.)**modifyxpw1) then
               retry = 0
            else
               retry = 1
            endif
         elseif(modifyx .eq. 6 ) then
!             drop small xxx
            xxx = xPimax* xNmax
            if(xxx .gt. 0.05) then
               retry = 0
            elseif( u .lt. (1.-xxx/0.05)**modifyxpw1 ) then
               retry = 0
            else
               retry = 1
            endif
         else
            retry = 0
         endif
      else
         retry = 0
      endif
      end


      subroutine chAcolx(pj, ia, iz, xs, a, ntp)
!          ????      
!      if only.. is omitted "both" in here and cheavyInt.f
!      at link time of large program, a message something like
!      TargetXs for x86-64 is missing
!      comes out and stops. 

!      use modColInfo
      implicit none
#include  "Zcondc.h"
#include  "Zcode.h"
#include  "Zptcl.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"


!
      type(ptcl):: pj   ! input  projectile ptcl
      integer ia    ! input. mass no. of target
      integer iz                ! input. charge no. of target
      real(8),intent(in):: xs ! x-sction on nucleus A=ia. mb
      type(ptcl):: a(*)  !  output. produced ptcls
      integer ntp   ! number of produced ptcls

#if defined (MacIFC)
!   if modColMediaInfo dose not work, activate next      
!  #include "Zworkaround.h"  
#endif


!
!c      if(IntModel .eq. 'int1') then
!c        if(pj.fm.p(4)  .lt. Elund) then
!            h-A collision  Fermi momentum is considered inside
      if( ActiveMdl .eq. 'phits') then
         call cphits(pj, ia, iz, xs, a, ntp)
      elseif( ActiveMdl .eq. 'dpmjet3') then
         call cdpmjet(pj, ia, iz, a, ntp)
      elseif( ActiveMdl .eq. 'jam') then
!                  TargetXs is the xs on 'ia'. This is used 
!            in jam to get bmax.  better to include it in 
!            chAcol argment in future ? now done so.
         call cjamEvent(pj, ia, iz, xs, a, ntp)
      elseif( ActiveMdl .eq. 'qgsjet2') then
         call cQGSjet(pj, ia, iz, a, ntp)
      elseif( ActiveMdl == 'epos') then
         call ceposGenOneEvent(pj, ia, iz, a, ntp)
      elseif( ActiveMdl == 'sibyll') then
         call csibyllevent(pj, ia, iz, a, ntp)
      elseif( ActiveMdl .eq. 'fritiof1.6' .or. 
     *        ActiveMdl .eq. 'nucrin') then
            call chALund(pj, ia, iz, a, ntp)
      elseif( ActiveMdl .eq. 'ad-hoc') then
         call chAcolAdhoc(pj, ia, iz, a, ntp)
!      elseif( ActiveMdl .eq. 'fritiof7.02') then
!         if(pj.code .eq. kkaon .and. LundPara(5) .ne. 0) then
!            call chAcolAdhoc(pj, ia, iz, a, ntp)
!         elseif(pj.code .eq. keta) then
!            call chAcolAdhoc(pj, ia, iz, a, ntp)
!         else
!            call chANewLund(pj, ia, iz, a, ntp)
!         endif
      elseif(ActiveMdl .eq. 'gheisha') then
!         if(pj.fm.p(4) .lt. Elund) then
         call chAGheisha(pj, ia, iz, a, ntp)
!         elseif(pj.fm.p(4)  .lt. Elund2) then
!            call chAcolAdhoc(pj, ia, iz, a, ntp)
!         elseif(pj.fm.p(4)  .lt. Elund3) then
!            if(pj.code .eq. kkaon .and. LundPara(5) .ne. 0) then
!               call chAcolAdhoc(pj, ia, iz, a, ntp)
!            elseif(pj.code .eq. keta) then
!               call chAcolAdhoc(pj, ia, iz, a, ntp)
!            else
!               call chANewLund(pj, ia, iz, a, ntp)
!            endif
      elseif(ActiveMdl .eq. 'incdpm3') then
!             dpmjet inclusive treatment
         call cincdpm3(pj, ia, iz, a, ntp)
      elseif(ActiveMdl .eq. 'byenergy') then
         if( pj%fm%p(4)-pj%mass .gt. 4.0) then
            call chAcolAdhoc(pj, ia, iz, a, ntp)
         else
            call chALund(pj, ia, iz, a, ntp)
         endif
      else
         call cerrorMsg(ActiveMdl, 1)
         call cerrorMsg('above model not supported',0)
      endif

      end
