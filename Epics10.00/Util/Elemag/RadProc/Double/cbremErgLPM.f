!     test cbremErgLPM
!
      real*8 ee, rho, eg
      integer i
      read(*,*) ee, rho
      ee = ee/1.e9      ! GeV
      rho = rho * 1.e3   ! kg/m3
      do i = 1, 10000
         call cbremErgLPM(ee, rho, eg)
         write(*, *) sngl(eg/ee)
      enddo
      end
      
!       sampling  Brems
!

      subroutine cbremErgLPM(ee, rho, eg)
      implicit none
      real*8 ee  !  electron energy in GeV
      real*8 rho !  air density in kg/m3
      real*8 eg  !  sampled gamma ray energy
!
      real*8  u, v, temp, s
      real*8 cmigdG, cmigdPsi, cmigdSforBrem, cmigdGzai

      real*8 a1, a2, vm, logvm, suma, a2prob
      parameter (vm = 1.d-4, logvm = -4*2.3025851, 
     *           a1 = 0.5, a2 = -4.d0*logvm/3.,
     *           suma = a1 + a2, a2prob = a2/suma)
!
      do while (.true.)
         call rndc(u)
         if(u .lt. a2prob) then
!             sampling by  1/v
            call rndc(u)
            v = exp(u * logvm)  ! fractional energy
!            rejection
            call rndc(u)
            if(u .lt. 1.-v) then
!                further check; get s
               s = cmigdSforBrem(ee, rho, v)
               if(u .lt. (1.-v) * cmigdGzai(s)*cmigdPsi(s)) then
!                  accepted
                  eg = v * ee
                  goto 100
               endif
            endif
         else
!              sampling by 2vdv
            call rndc(u)
            call rndc(v)
            v = max(u, v)
!            rejection, get s
            s = cmigdSforBrem(ee, rho, v)
            temp = cmigdGzai(s) * (cmigdG(s) + 2* cmigdPsi(s))/3.
            call rndc(u)
            if(u .lt. temp) then
!               accepted
               eg = v * ee
               goto 100
            endif
         endif
      enddo
 100  continue
      end
