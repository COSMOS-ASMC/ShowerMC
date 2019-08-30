      subroutine cfixProc(id, rl2kg)
!     select shortest one  among sampled interaction lengths.
!     First select shortest one among those interactions of the
!     length is given in kg/m2 or r.l (not in actual length)
!     then convert it to the atcual length in m.
!     Then compare it with other lengths in m and take shortest one.
!     We do this because to convert kg/m2 or r.l into m needs some numerical
!     compuation time  due to gradual change of density.

      use modSetIntInf
      use modAtmosDef
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackv.h"
! #include  "Zearth.h"

      integer,intent(in):: id   ! 1,2 --> projectile is g / e
                                ! >= 3 --> //                mu, hadron

      real(8),intent(in):: rl2kg ! r.l x rl2kg --> kg/m2.
            ! this is media%X0kg  to be used when id =1 or 2.

      real*8 h, leng, t, minlen, pcut, tk
      real(8),external:: clen2thick
!   
      integer  jcut

      integer:: idxt, idxm
      real(8):: mint, minm

!     find min one among those with "thickness"
!     (i.e, exclude  decay, mag{pair,brem})
      ProcessNo =  0
      call cfixProc0(idxt, mint)    ! output mint is thickess in  kg/m2 or  r.l
                              ! EM--> r.l  Hadint -> kg/m2  
      call cfixProc1(idxm, minm)     ! those given in m
      
      h = TrackBefMove%pos%height ! used later
      if( idxt  >  0 ) then
         if( mint < Infty ) then
          !   convert it to m ; first make it kg/m2 --> tk
            if( id > 2 ) then
               tk =   mint      ! in kg/m2
            else
               tk =   mint * rl2kg ! r.l --> kg/m2
               IntInfArray(idxt)%thickness = tk ! kg/m2 value needed
                                       ! for Energy loss calc.
            endif
         else
            tk = Infty
         endif
!         kg/m2 to length in m.
         call cthick2len(TrackBefMove, tk, minlen, t, jcut)
         if( idxm > 0 ) then
            if( minm <= minlen ) then
               ProcessNo = idxm
               IntInfArray(idxm)%thickness = clen2thick(h,
     *              TrackBefMove%vec%coszenith,
     *              IntInfArray(idxm)%length )
               MoveStat = ToInteract
            else
               ProcessNo = idxt
            endif
         else
            ProcessNo = idxt
         endif
         if( ProcessNo == idxt ) then
            if( jcut /=0 ) then
!               too long tk. so cut to t
!             IntInfArray(idxt)%process = ' '  ! keep the process name
!                    though the no such process will
!                    take place by the next MoveStat.
               
               MoveStat = Truncated
               IntInfArray(idxt)%thickness = t ! its length is minlen
            endif
            IntInfArray(idxt)%length = minlen !
         endif
      elseif( idxm > 0 ) then
!            idxt < 0
         ProcessNo = idxm
         IntInfArray(idxm)%thickness = clen2thick(h,
     *              TrackBefMove%vec%coszenith,
     *              IntInfArray(idxm)%length )
      endif
      if( ProcessNo == 0 ) then
         ProcessNo = 1
         IntInfArray(ProcessNo)%length = Infty
!         IntInfArray(ProcessNo)%thickness = Infty         
      endif
!         In the case  of muon, if individual knockon process  
!       is neglected (by parameter setting or with high Emin)
!       the length could be quite large (say, 6000  km).
!       and results in error in the next call.
!       To avoid that, we cut the path here
!
      if( TrackBefMove%vec%coszenith .lt. 0.) then
!     pcut = 300.d3
         pcut = Eradius / 20d0
      else
!     pcut = 30.d3
         pcut = Eradius / 200d0
      endif
      if(IntInfArray(ProcessNo)%length .gt. pcut) then
         MoveStat = Truncated
         IntInfArray(ProcessNo)%length  = pcut
         IntInfArray(ProcessNo)%thickness = clen2thick(h, 
     *      TrackBefMove%vec%coszenith, 
     *      IntInfArray(ProcessNo)%length )
      endif
      end
      subroutine cfixProc0(idxt, mint)
!         select min thickness among non "decay" processes
      use modSetIntInf
      implicit none
#include  "Zglobalc.h"
      integer,intent(out):: idxt ! 0 if no data. else idxt-th
              ! process is minimum thikness
      real(8),intent(out):: mint
      integer:: i

      idxt = 0
      mint = Infty 
!          simply check non "decay" (thickness in r.l)  part
      do i = 1, NumberOfInte
         if(.not. IntInfArray(i)%decay) then
            if( IntInfArray(i)%thickness < mint ) then
               mint = IntInfArray(i)%thickness
               idxt= i
            endif
         endif
      enddo
      end

      subroutine cfixProc1(idxm, minm)
!     select min length among "decay" processes
!     (decay=T means   IntInfArray(i)%length has been
!     given in m. Process may be not only decay; it
!     may be one of magnetic brems, magnetic pair reation..
!      
      use modSetIntInf
      implicit none
#include  "Zglobalc.h"
      integer,intent(out):: idxm
      real(8),intent(out):: minm
      integer:: i

      idxm = 0
      minm = Infty 
!          simply check "decay" (actual length)  part
      do i = 1, NumberOfInte
         if(IntInfArray(i)%decay) then
            if( IntInfArray(i)%length < minm ) then
               minm = IntInfArray(i)%length
               idxm = i
            endif
         endif
      enddo
      
      end
