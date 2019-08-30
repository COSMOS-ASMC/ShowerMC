!     test epGetXXsec
      implicit none
#include "Zmedia.h"
#include "ZepManager.h"
       type(epmedia)::   media
      integer  icon, n, i
!
      character*8 name
      real*4  xsec(7), Ex
      MediaDir(1)="$EPICSTOP/Data/Media/"
      write(0,*) "Enter media name"
      read(*,'(a)')  media%name 
      write(0,*) "media name=", media%name
      call epReadXXsec(media, icon)
      write(0,*) ' n=', media%xcom%size, ' icon =',icon
!

      do i = 1, media%xcom%size
         write(0,'(1p,7g14.4)') media%xcom%tab(:,i)
      enddo
!
!     subroutine cGetXXsec(Ex, xcomtab, n, m1, m2, xsec, icon)
!      implicit none
!        give X-ray xsection;  1) coherent scatt
!                              2) incoherent (compton) scatt.
!                              3) photo-absorption,
!                              4)  pair prod.by nucl
!                              5)  pair prod. by atomic elec.
!                              6)  attn coef. with coh.
!                              7)  attn coef. without coh.
      Ex = 1.e-6
      do i = 1,700
         call cGetXXsec(Ex, media%xcom%tab,
!     *   media.xcom.size, 7, xsec, icon)
     *   media%xcom%size, 1, 7, xsec, icon)
         if( icon /= 0 ) then
            write(0,*) ' cGetXXsec icon =',icon
            stop
         endif

         write(*,'(1p,8g14.4)') Ex, xsec(1:7)
         Ex = Ex*10.**0.01
      enddo
      end
