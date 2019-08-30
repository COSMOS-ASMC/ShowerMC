!     test epGetXXsec
      implicit none
      integer  icon, n

      character*8 name
      real*8 xcomtab(8, 200)
      name = 'BGO'
      call epReadXXsec(name, xcomtab, n, icon)
      write(0,*) ' n=', n, ' icon =',icon
      end

      subroutine epGetXXsec(Ex, xcomtab, n, m, xsec, icon)
      implicit none
!        give X-ray xsection;  1) coherent scatt
!                              2) incoherent (compton) scatt.
!                              3) photo-absorption,
!                              4)  pair prod.by nucl
!                              5)  pair prod. by atomic elec.
!                              6)  attn coef. with coh.
!                              7)  attn coef. without coh.
      real*8  Ex  ! input.  X-ray energy in GeV. 1keV to 100GeV
      integer   n    ! input. row size of xcomtab. the column is 8
      real*8  xcomtab(8, n)  ! input.  x-section table obtained by using 
                     ! CreateTab   and stored as bgo.xcom etc
                     ! values are in log. (original 0 --> -100)
                     ! (1,n)=E(GeV), (2,n)=coh. (3,n)=incoh. (4,n)=p.e.
                     ! (5,n)=n.pair (6,n)=e.pair, (7,n)=atten(with coh)
                     ! (8,n)=atten.(without  coh)
      integer   m    ! input. first m xsections are obtained in xsec. in
                     !        unit of 1/(g/cm^2)
      real*8   xsec(8)  ! output.  at least size m.  xsec(k)  is k-th xsection.
      integer icon   ! output.  0--> ok.  1-->Ex<1keV. extrapolation not guaranteed
                     !          2--> Ex>100 GeV.  values at 100 GeV is given
      real*8   Exl
      real*8   dx, grad
      integer i, loc

      Exl = log(Ex)
      if( Exl .lt. xcomtab(1, 1) ) then
         do i = 1, m
            xsec(i) = xcomtab(i+1, 1)
         enddo
      elseif(EXl .gt. xcomtab(1,n) ) then
         Exl = xcomtab(1,n)
      else
!          find i(=loc) such that   Ei <= Eg < Ei+1  ( i=1, 2, ...n-1)
         call kdwhereis(Exl, n, xcomtab, 8, loc)
!           if Ex=100, loc = n
         if( loc .lt. n ) then
            if( xcomtab(1, loc) .eq.  xcomtab(1, loc+1) ) then
               loc = loc + 1
            endif
            if(loc .eq. n) then
               do i = 1, m
                  xsec(i) = xcomtab(i+1, n)
               enddo
            else
               dx = xcomtab(1, loc+1)- xcomtab(1, loc) 
               if(dx .eq. 0.) then
                  do i = 1, m
                     xsec(i) = xcomtab(i+1, loc)
                  enddo
               else
                  do i = 1, m
                     grad =(xcomtab(i+1, loc+1)-xcomtab(i+1, loc))/dx 
                     if(i .eq. 4 .or. i .eq. 5) then
                        if(xcomtab(i+1, loc+1) .eq. -100.) then
                           xsec(i)= -100.
                        else
                           xsec(i) = grad* (Exl- xcomtab(1, loc)) 
     *                       + xcomtab(i+1, loc)
                        endif
                     else
                        xsec(i) =  grad * (Exl- xcomtab(1, loc)) 
     *                          + xcomtab(i+1, loc)
                     endif
                  enddo
               endif
            endif   
         else
            do i = 1, m
               xsec(i) = xcomtab(i+1, n)
            enddo
         endif
      endif
      do i = 1, m
         if(xsec(i) .eq. -100) then
            xsec(i)  = 0.
         else
            xsec(i) = exp(xsec(i))
         endif
      enddo
      end

      subroutine epReadXXsec(name, xcomtab, n, icon)
      implicit none
#include "ZepManager.h"
      character*8 name   ! input  media name.
      real*8  xcomtab(8,*)  ! output.  E vs xsection  tab
      integer n             ! output.  number of rows of xcomtab
      integer icon          ! output.  icon = 0.  xcombtab obtained
                            !  -1.  no. tab found.

      integer klena
      integer i, j, result, k
      character*150 mdpath
      
      icon = -1
      do i = 1, MaxMediaDir
         if( klena(MediaDir(i)) .gt. 0 ) then
            mdpath = ' '
            mdpath = MediaDir(i)(1:klena(MediaDir(i)))
     *            // '/'  //name(1:klena(name))//".xcom"
            call copenf(iowk, mdpath, result)
            if( result .ne. 0) then
               write(0,*) mdpath
               write(0,*) ' cannot be opended'
               stop 999
            endif
            call cskipComment(iowk, icon)
            if( icon .ne. 0) then
               write(0,*)
     *              "#------------- line not found in .xcom"
               stop 987
            endif
            n = 0
            do while(.true.)
               read(iowk, *, end=100) ( xcomtab(k,n+1), k = 1, 8)
               n = n + 1
            enddo
 100        continue
            icon = 0
            do k = 1, n
               do j= 1, 8
                  if( xcomtab(j, k) .eq. 0.) then
                     xcomtab(j,k) = -100
                  else
                     xcomtab(j,k) = log( xcomtab(j,k) )
                  endif
               enddo
            enddo
            exit
         endif
      enddo                               
      end
