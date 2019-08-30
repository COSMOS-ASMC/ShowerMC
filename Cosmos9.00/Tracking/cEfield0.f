!     defninition of modEfield is in LibLoft/Module/cepControl.f
!     Next two sub are used only at test time
!        
      subroutine cEfieldRead
      use modEfield
      implicit none
      read(*, Efparam) 
      end subroutine cEfieldRead

      subroutine cEfieldWrite
      use modEfield      
      implicit none
        write(0, Efparam) 
      end      subroutine cEfieldWrite

      subroutine cEfieldCheck
      use modEfield      
      implicit none
      integer:: i

      if(HowEfield /= 1 ) return  !!!!!
      do i =  nmaxEfield, 1, -1   ! last non zero Ef 
         if(dot_product(myEf(i)%Ef(:),myEf(i)%Ef(:))  > 0 ) then
            nEfield = i
            exit
         endif
      enddo
      checkH(1:nEfield) =
     *     myEf(1:nEfield)%H1**2+myEf(1:nEfield)%H2**2 > 0.
      checkR(1:nEfield) =
     *     myEf(1:nEfield)%R1**2+myEf(1:nEfield)%R2**2 > 0.
      checkT(1:nEfield) =
     *     myEf(1:nEfield)%T1**2+myEf(1:nEfield)%T2**2 > 0.
      do i = 1, nEfield
         useH = useH .or. checkH(i)
         useR = useR .or. checkR(i)
         useT = useT .or. checkT(i)
      enddo
      useH = useH .and. .not. useT  ! useT has priority
      
      if( nEfield >= 1 ) then
              ! check if height is from low to heigh
         do i = 1, nEfield
            if( myEf(i)%H1 > myEf(i)%H2) then
               write(0,*)  i,'-th Efiled height error'
               write(0,*)
     *              'H1=', myEf(i)%H1,' >', myEf(i)%H2,'=H2'
               stop
            endif
            if( myEf(i)%R1 > myEf(i)%R2) then
               write(0,*)  i,'-th Efiled Radius error'
               write(0,*)
     *              'R1=', myEf(i)%R1, ' >',  myEf(i)%R2,'=R2'
               stop
            endif
            if( myEf(i)%T1 > myEf(i)%T2) then
               write(0,*)  i,'-th Efiled Time error'
               write(0,*)
     *              'T1=', myEf(i)%T1, ' >',  myEf(i)%T2,'=T2'
               stop
            endif
            if( nEfield > i ) then
               if( myEf(i+1)%H1 < myEf(i)%H2 ) then
                  write(0,*)
     *                 i,'-th Ef height region<=',myEf(i)%H2
                  write(0,*)
     *                 ' but is not < next region=',myEf(i+1)%H1
                  stop
               endif
            endif                  
         enddo
      endif
      write(0,*) ' # of Efield(s)=',nEfield
      end subroutine cEfieldCheck


      subroutine cEfieldDet2xyz
      ! convert Ef in det system into E-xyz system
      ! should be called once before starting tracking
      use modEfield
      implicit none
      integer::i
      do i = 1, nEfield
         call cdet2xyzD(myEf(i)%Ef(:), Efxyz(:,i))
      enddo
      end  subroutine cEfieldDet2xyz

      subroutine cEfield0(aTrack, Efout)
!      Electric field is given depending on particle position
!      and/or time.  
!      H          T                      R
!    given   not given --> H is used   R at that H is used. if R is out 
!                                      of range, E=0. else E at H
!                                      if R is not given  E at H
!  not given    given  --> T is used   same as above, use T instead of H.
!    given      given  --> T is used    //
!  not gvien   not given               use first R,.if it exists
!                                      if in range, use E(1) else E=0
!     if R absent, use E(1)
      use modEfield        
      implicit none
#include  "Zglobalc.h"
#include "Ztrack.h"
! #include "Zmagfield.h"
#include "Zincidentv.h"
      type (track):: aTrack  ! input. Track info in E-xyz.
         ! position and "where" info are used.
      real(8),intent(out):: Efout(3) ! Elecric field
                       ! in E-xyz system.  V/m.
      integer:: i, hindex, tindex, out
      logical,save:: first=.true.
      real(8)::time
      
      real(8)::r
      real(8)::Zero(3)=(/0.,0.,0./)

      if( first ) then   ! convert Ef in det into E-xyz
         call cEfieldCheck
         call cEfieldDet2xyz
         first = .false.
      endif

      hindex = 0
      if(useH) then
         do i = 1, nEfield
            if(checkH(i)) then
               if( aTrack%pos%height >= myEf(i)%H1 .and.
     *              aTrack%pos%height <= myEf(i)%H2 ) then
                  hindex = i
                  exit
               endif
            endif
         enddo
         if( hindex == 0) then  ! H is out of range
            Efout(:) = Zero(:)
         else
              !  see if R is in the range
            call cEfieldCheckR(aTrack, hindex, out)
            if(out == 0 ) then
             !  Efout(:) = myEf(hindex)%Ef(:)  ! bug
               Efout(:) = Efxyz(:,hindex)
            else
               Efout(:) = Zero(:)
            endif
         endif
      elseif(useT) then
         tindex = 0
         do i = 1, nEfield
            if(checkT(i)) then
               time = aTrack%t/c*Tonsec
               if( time >= myEf(i)%T1 .and.
     *             time <= myEf(i)%T2 ) then
                  tindex = i
                  exit
               endif
            endif
         enddo
         if( tindex == 0 ) then  ! T is out of range
            Efout(:) = Zero(:)
         else
              !  see if R is in the range
            call cEfieldCheckR(aTrack, tindex, out)
            if(out == 0 ) then
!               Efout(:) = myEf(tindex)%Ef(:)  bug
               Efout(:) = Efxyz(:,tindex)
            else
               Efout(:) = Zero(:)
            endif
         endif
      else  ! H,T is not used. 
         i = 1   ! use first data if possible
         if( checkR(i) ) then   ! check first R if any
            call cEfieldCheckR(aTrack, i, out)
            if(out == 0 ) then
!               Efout(:) = myEf(i)%Ef(:)
               Efout(:) = Efxyz(:,i)
            else
               Efout(:) = Zero(:)
            endif
         else  ! no check of R.  assume first Ef is applied
!            Efout(:) = myEf(i)%Ef(:)
            Efout(:) = Efxyz(:,i)
         endif
      endif
      end   subroutine cEfield0
      
      subroutine cEfieldCheckR(aTrack, i, out)
      use modEfield
      implicit none
#include "Ztrack.h"
!  #include "Zmagfield.h"
      
      type (track):: aTrack  ! input
      integer,intent(in):: i ! R or H index
      integer,intent(out)::out ! 0 aTrack is in the range or no check is needed
                   !  else out of range
      real(8):: r

      if( checkR(i) ) then
                ! core distance
         if(DefofR == 'p' ) then ! perpedicular distance
            call cgetPcoreDist(aTrack, r)
         elseif(DefofR == 'h') then ! horizontal distance
            call cgetHcoreDist(aTrack, r)
         else
            write(0,*) ' DefofR =',DefofR, ' invalid'
            stop
         endif
         if( r >= myEf(i)%R1 .and.
     *        r <= myEf(i)%R2 ) then
            out = 0
         else                   ! R is out of range
            out = 1
         endif
      else
         out = 0
      endif
      end subroutine cEfieldCheckR

      subroutine cgetPcoreDist(aTrack, r)
      implicit none
#include "Ztrack.h"
! #include "Zmagfield.h"      
#include "Zincidentv.h"
#include  "Zobs.h"
#include "Zobsv.h"
      type (track):: aTrack  ! input. in E-xyz system
      real(8),intent(out):: r ! Perpendicular core distance

      real(8):: rdet(3)  !   ptcl pos in det-sys

      call cxyz2det(ObsSites(aTrack%where)%pos%xyz,
     *              aTrack%pos%xyz%r(:), rdet(:))
      r = abs(dot_product(rdet(:), AngleAtObsCopy%r(:)))  ! Angle..  is dwn.g
      end      subroutine cgetPcoreDist

      subroutine cgetHcoreDist(aTrack, r)
      implicit none
#include "Ztrack.h"
! #include "Zmagfield.h"
#include "Zincidentv.h"
#include "Zobs.h"
#include "Zobsv.h"
      type (track):: aTrack  ! input. in E-xyz system
      real(8),intent(out):: r ! Perpendicular core distance

      real(8):: rdet(3)  !   ptcl pos in det-sys
      real(8):: Raxis(3) ! shower core at the same height
             !         as the ptcl


      call cxyz2det(ObsSites(aTrack%where)%pos%xyz,
     *               aTrack%pos%xyz%r(:), rdet(:))
      Raxis(:) = 
     *  (aTrack%pos%height - ObsSites(aTrack%where)%pos%height) /
     *  AngleAtObsCopy%r(3)     ! length from origin( <0 )
     *  * AngleAtObsCopy%r(:)   ! vector (now upgoing) since
                       ! AngleAtObsCopy is downogin.
      r = sqrt(sum( ( rdet(:)-Raxis(:))**2 ))
      end      subroutine cgetHcoreDist

