!     *****************************
      subroutine epDraw_wave(comp, p, n)
      implicit none
#include "Zglobalc.h"
#include "Zep3Vec.h"
#include "Zcnfig.h"
#include "ZepDraw.h"
       type(Component)::  comp  ! input. component
       type(epPos)::  p(*)     ! output. (x,y,z) to describe
                              !   an wave in local coord.
                              ! (x,y,z)= gpsep is a separator
                              ! to be converted to a blank line

      integer  n              ! output.  number of (x,y,z) data
                              ! put in p.



      logical kdgtest, pused

      integer nvhccl
      parameter (nvhccl = nvccl/2+4)
       type(epPos)::  ptemp(nvhccl)

      real*8  x, yat, cenx, temp

      integer m, i, nw, j, j1, j2, j3
      integer  nmin, nmax      
      real*8 t1, t2
      real*8  eps
      save eps


      real*8 ep_wavef
      external ep_wavef
!     *************************************
      real*8 L, h, d, s, e, a, hL, qL, t
      common /Zwave/  L, h, d, s, e, a, hL, qL, t
!     *************************************

      data eps/1.d-9/
!         center of the m-th half wave
      cenx(m) = m*hL + qL
!         y value at x of the circle centered at cenx(m)
      yat(x, m) = (-1)**m * sqrt( abs( qL**2 - (x-cenx(m))**2 ))
!
!         set common area      
      call ep_wavesetcom(comp)

      n =0
      nw = 0
      nmin =  int( (s+eps)/hL )
      nmax =  int( (e-eps)/hL )

      if( kdgtest(how, 1) .or. kdgtest(how, 2) .or.
     *    kdgtest(how, 5) .or. kdgtest(how, 6) ) then
!         we first compute upper-front wave coord.
!         and save them in p(1) ~ p(nw) (including last
!         sepaator).  These may be modified later, if
!         only kdgest(how,1) is on.
!          upper wave 
         do i = nmin, nmax
            if( i  .eq. i/2*2 ) then
               if( i .eq. nmin ) then
                  t2 =
     *            atan2( yat(s+eps, nmin), s+eps -cenx(nmin) )
     *            *Todeg
                  if(t2 .le. 0.) then
                     t2 = t2 + 360.d0
                  endif
               else
                  t2 = 180.d0
               endif
               if( i .eq. nmax ) then
                  t1 = 
     *            atan2( yat(e-eps, nmax), e-eps-cenx(nmax))
     *            *Todeg
                  if(t1 .le. 0.d0) then
                     t1 = -t1
                  endif
               else
                  t1 = 0.d0
               endif
            else
               if( i .eq. nmin ) then
                  t1 = atan2( yat(s, nmin), s-cenx(nmin))*Todeg
                  if(t1 .lt. 0.d0) then
                     t1 = 360.d0 + t1
                  endif
               else
                  t1 = 180.d0
               endif
               if( i .eq. nmax) then
                  t2 = atan2( yat(e, nmax), e-cenx(nmax))*Todeg
                  t2 = min(360.d0 + t2, 360.d0)
               else
                  t2 = 360.d0
               endif
            endif

            if(t1 .gt. t2 ) then
               temp = t1
               t1 = t2
               t2 = temp
            endif

            call epdrawCcl(qL, 0.d0, t1, t2, ptemp, m)


            if(i .eq. nmin) then
               if(i .eq. i/2*2) then
                  j1 = m-1
                  j2 = 1
                  j3 = -1
               else
                  j1 = 1
                  j2 = m-1
                  j3 = 1
               endif
            else
               if(i .eq. i/2*2) then
                  j1 = m-1
                  j2 = 2
                  j3 = -1
               else
                  j1 = 2
                  j2 = m-1
                  j3 =1
               endif
            endif
            do j = j1, j2, j3
               n = n + 1
               p(n)%x = ptemp(j)%x + cenx(i)
               p(n)%y = 0.
               p(n)%z = ep_wavef(comp, 1,  p(n)%x )
            enddo
         enddo

         n = n +1
         p(n)%x = gpsep

         nw = n
      endif

      if( kdgtest(how,1) .and. (.not. kdgtest(how,2) .and.
     * .not. kdgtest(how,5) .and. .not. kdgtest(how,6))) then

!         only bottom wave is specified
!         change the p(1),p(nw)
 
        do i = 1, nw - 1
           p(i)%z = p(i)%z - t
        enddo
         
         
        do i = 1, nw - 1
           n = n +1 
           p(n)%x = p(i)%x
           p(n)%y = p(i)%y + d
           p(n)%z = p(i)%z
        enddo
        n = n + 1
        p(n)%x = gpsep
        
        n = n + 1
        p(n)%x = gpsep
        goto 100
      endif

      if( kdgtest(how, 5) .and. (.not. kdgtest(how,1) .and.
     *   .not. kdgtest(how,2) .and. .not. kdgtest(how,6))) then
!          only 5 is specified.  p must be modified
         do i = 1, nw - 1
            p(i)%y = d
         enddo
         do i = 1, nw - 1
            n = n + 1
            p(n)%x = p(i)%x
            p(n)%y = p(i)%y
            p(n)%z = p(i)%z - t
         enddo
         n = n + 1
         p(n)%x = gpsep
        
         n = n + 1
         p(n)%x = gpsep
         goto 100
      endif

      if( kdgtest(how, 6) ) then
!           upper wave
!           at y = d.   p(1) ~ p(nw) should be front coord. already
         do i = 1, nw -1
            n = n +1
            p(n)%x = p(i)%x
            p(n)%y = d
            p(n)%z = p(i)%z
         enddo
         n = n + 1
         p(n)%x = gpsep

         n = n + 1
         p(n)%x = gpsep
         pused = .true.
!        =============
      endif

      if( kdgtest(how, 1)) then
!     
         if(pused) then
            do i = 1, nw-1
               n = n + 1
               p(n)%x = p(i)%x
               p(n)%y = p(i)%y + d
               p(n)%z = p(i)%z - t
            enddo

            n = n + 1
            p(n)%x = gpsep

            do i = 1, nw-1
               n = n + 1
               p(n)%x = p(i)%x
               p(n)%y = p(i)%y
               p(n)%z = p(i)%z - t
            enddo

         else
!             we should use p(1)...
            do i = 1, nw - 1
               p(n)%x= p(i)%x
               p(n)%y = d
               p(n)%z = p(i)%z
            enddo
            pused = .true.
         endif
         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
      endif
!        =============

      if( kdgtest(how, 2) ) then
         if(pused) then
            do  i = 1, nw -1
               n = n + 1
               p(n)%x = p(i)%x
               p(n)%y = p(i)%y
               p(n)%z = p(i)%z - t
            enddo
            n = n + 1
            p(n)%x = gpsep

            do i = 1, nw 
               n = n + 1
               p(n) = p(i)
            enddo
         else
!          we should used p(1)..
            do i = 1, nw -1
               n = n+1
               p(n)%x = p(i)%x
               p(n)%y = 0.
               p(n)%z = p(i)%z - t
            enddo
            pused = .true.
            n = n + 1
            p(n)%x = gpsep
         endif
         n = n + 1
         p(n)%x = gpsep
      endif
         
      if( kdgtest(how, 5) ) then
!         if(pused) then  ! should be always so here
         do i = 1, nw -1
            n = n + 1
            p(n)%x = p(i)%x
            p(n)%y = p(i)%y + d
            p(n)%z = p(i)%z
         enddo
         n = n + 1
         p(n)%x = gpsep
         
         do  i = 1, nw -1
            n = n + 1
            p(n)%x = p(i)%x
            p(n)%y = p(i)%y + d
            p(n)%z = p(i)%z - t
         enddo
!         else
!            should be error            
!         endif
         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
      endif

 100  continue
      if( kdgtest(how, 3) )then
         n = n + 1
         p(n)%x = s
         p(n)%y = 0.
         p(n)%z = ep_wavef(comp, 1, s)
         n = n +1
         p(n)%x = s
         p(n)%y = d
         p(n)%z = p(n-1)%z 
         n = n + 1
         p(n)%x = gpsep

         n = n + 1
         p(n)%x = s
         p(n)%y = 0.
         p(n)%z = ep_wavef(comp, 0, s)

         n = n +1
         p(n)%x = s
         p(n)%y = d
         p(n)%z = p(n-1)%z 

         n = n + 1
         p(n)%x = gpsep

         n = n + 1
         p(n)%x = gpsep
!        =============
      endif

      if( kdgtest(how, 4) )then
         n = n + 1
         p(n)%x = e
         p(n)%y = 0.
         p(n)%z = ep_wavef(comp, 0, e)

         n = n +1
         p(n)%x = e
         p(n)%y = d
         p(n)%z = p(n-1)%z 

         n = n + 1
         p(n)%x = gpsep


         n = n + 1
         p(n)%x = e
         p(n)%y = 0.
         p(n)%z = ep_wavef(comp, 1, e)

         n = n +1
         p(n)%x = e
         p(n)%y = d
         p(n)%z = p(n-1)%z 
         n = n + 1
         p(n)%x = gpsep
         n = n + 1
         p(n)%x = gpsep
      endif
      end
