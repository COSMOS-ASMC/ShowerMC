!c     *****************
!      subroutine eptotcbS(vmin, vmax, ans)
!c     ****************
!*     implicit none
!
!      real*8 vmin, vmax, ans
!c
!c     vmin and vmax  are Eg/Ee.
!c     ans is in mb
!c
!c        integration of Seltzer's bremsung function from vmin to vmax.
!c        vmin > 0. 
!c        see epBrSfs for preparation.
!c
!      external  epBrSfs
!      real*8    epBrSfs
!
!      real*8 v1, v2,  ans1
!
!      ans=0.
!      v2=vmax
!c      v1 = 0.
!c////////////
!      v1 = max( v2*0.9d0,  vmin)
!      call k16pGaussLeg(epBrSfs, v1, v2, 16,  ans)
!      v2 = v1
!      do while (v1 .ne. vmin)
!         v1 = max( v2/2.d0,  vmin)
!         call k16pGaussLeg(epBrSfs, v1, v2, 16,  ans1)
!         ans=ans+ans1
!         v2=v1
!      enddo
!      end
!
!      **********************************
      function epBrSfs(media, Eeme, x) result(ans)
!
!     function value:  cross-section  ds/dx in mb
!      x = Eg/Ee
!         Before calling this, common /BPgene/ must be
!       fixed, and epZforSeltz must be called.
!
!       
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
       type(epmedia)::  media    ! input
      real(8),intent(in):: Eeme  ! Ee/me
      real(8),intent(in)::x   ! Eg/Ee
      real(8)::ans  !  ds/dx in mb
!
!   
      integer i
      real*8 sum, Ek, beta2, Ee 
      real*8 epBrSf, k

      Ee = Eeme*masele
      Ek = Ee - masele
      k = Ee*x/Ek
      beta2  =  1.d0 -  1.d0/(Ek/masele +1)**2
      sum =0.
      do i = 1,  media%noOfElem
!         epBrSf= (beta/Z)**2 k ds/dk  ( in mb )
!         k=Eg/Ek ; x=Eg/Ee; k=x Ee/Ek; dk=dx Ee/Ek
!         1/dx= Ee/Ek/dk
!         ds/dx = ds/dk* Ee/Ek = (Z/beta)**2*epBrSf/k *Ee/Ek

         sum = sum + 
     *    (media%elem(i)%Z**2/beta2) 
     *         *  epBrSf(media%elem(i)%Z, Ek, k) /k *Ee/Ek 
     *         *  media%No(i)
      enddo
!            ds/dx has been obtained;
      ans = sum
      end function epBrSfs
!     ********************************************
      module SeltzerTab
      integer,parameter:: maxnoofz=30  ! max memorizable z vlaues
      integer,parameter:: vm = 30
      integer,parameter:: em = 57
      real(8),save:: dsdk(vm, em, maxnoofz)
      integer,save:: zvalue(maxnoofz)
      integer,save:: cnter
      end module SeltzerTab

      function epBrSf(Z, Ek, k) result(ans)
      use SeltzerTab
      implicit none
#include "Zmass.h"
      real*8 Z ! input Z
      real*8 Ek ! input  Electron kinetic energy in GeV
      real*8 k  ! input Eg/Ek
      real(8)::ans ! 
!       function value. (beta/Z)**2 k ds/dk in mb
!           x= Eg/Ek.
! NOTE:  Seltzer's table gives  (beta/Z)**2 *k ds/dk in mb 
!        So if you want to get ds/dx for x=Eg/Ee,
!        you must do following for the ans of this function
!        ds/dk =  (Z/beta)**2 ans/k
!        since kEk =xEe(=Eg), i.e, Ekdk=Eedx,  1/dx =Ee/Ek *  1/dk,
!        we have  ds/dx =  (Ee/Ek)* ds/dk = (Ee/Ek)* (Z/beta)**2 
!         *ans/k
! 

      integer  i
      integer intz 
      real*8  zintegral, f1, f2
      real*8  epBrSeltzer, Ee

      intz = Z
      zintegral = intz
      Ee = Ek + masele

      do i = 1, cnter
         if(zvalue(i) .eq. intz) then
            f1= epBrSeltzer(intz, dsdk(1,1,i), vm, em, Ek, k)
            goto 20
         endif
      enddo
      call cerrorMsg('no Z found in epBrSf',0)
 20   continue
      if(zintegral .ne. Z) then
         do i = 1, cnter
            if(zvalue(i) .eq. intz+1) then
               f2 = epBrSeltzer(intz+1, dsdk(1,1,i), vm, em, Ek, k)
               goto 30
            endif
         enddo
         call cerrorMsg('no Z found in epBrSf',0)
 30      continue
         ans = (f2-f1)*(Z - intz) + f1
      else
         ans = f1
      endif
!       function value. (beta/Z)**2  x ds/dx in mb 
!        x= Eg/Ek.
      end function epBrSf
!     ******************************
      subroutine epZforSeltz(ini, z)
      use SeltzerTab
      implicit none
!
!        inform z value to be used. This must be called before
!        epBrSf is used.
!
      integer ini  ! if 1, internal counter is cleared.
                   !  else, internal counter is increased and
                   !  z value is memorized together with previous
                   !  values. For the first call must be 0
      real*8 z     ! z value.  If not integral value, integer z1, z2 
                   !           such that z1 < z < z2 are memorized.
      
!
      real*8  zintegral
      integer intz


      if(ini .eq. 1) then
         cnter = 0
      endif
      intz = z
      zintegral = intz
      call epZcheck(intz)
      call epSeltzXsec(intz, dsdk(1,1,cnter), vm, em)
      if(z .ne. zintegral) then
         call epZcheck(intz+1)
         call epSeltzXsec(intz+1, dsdk(1,1,cnter), vm, em)
      endif
      end subroutine epZforSeltz
!     ***************
      subroutine epZcheck(intz)
      use SeltzerTab
      implicit none
      integer intz ! input Z

      integer i
      do i = 1, cnter
         if(intz .eq. zvalue(i)) goto 100
      enddo
      cnter = cnter + 1
      if(cnter .gt. maxnoofz) then
         call cerrorMsg('too many z to  epZforSeltz', 0)
      endif
      zvalue(cnter) = intz
 100  continue
      end subroutine epZcheck
!     **************************************************
      function
     *       epBrSeltzer(intz, dsdks, kmax, emax, Ekin, k)
     *       result( ans )
      use SeltzerTab
      implicit none
      integer intz   ! input  Z  value
      integer kmax, emax  ! input. see below
      real*8 dsdks(kmax, emax)  ! input cross-section table  kds/dk mb
      real*8 Ekin !  input. Electron kinetic energy in GeV.
      real*8 k  !  input. Eg/(Ee-me)
      real(8):: ans !  (beta/Z)** kds/dk in mb
!

      integer,parameter::points =4


      real*8  Ek(em), kv(vm)
   !       even array, Gfortran treats this automatic so
   !       we must use "save"  
      real(8),save:: Eklog(em)
      real*8  error
      logical,save::first= .true.

      integer:: i

        data Ek/
     * 1d-06,    1.5d-06,   2d-06,     3d-06,    4d-06,
     * 5d-06,    6d-06,     8d-06,     1d-05,    1.5d-05,
     * 2d-05,    3d-05,     4d-05,     5d-05,    6d-05,
     * 8d-05,    0.0001d0,  0.00015d0, 0.0002d0, 0.0003d0,
     * 0.0004d0, 0.0005d0,  0.0006d0,  0.0008d0, 0.001d0, 
     * 0.0015d0, 0.002d0,   0.003d0,   0.004d0,  0.005d0,
     * 0.006d0,  0.008d0,   0.01d0,    0.015d0,  0.02d0,
     * 0.03d0,   0.04d0,    0.05d0,    0.06d0,   0.08d0,
     * 0.1d0,    0.15d0,    0.2d0,     0.3d0,    0.4d0,
     * 0.5d0,    0.6d0,     0.8d0,     1d0,      1.5d0, 
     * 2d0,      3d0,       4d0,       5d0,      6d0, 
     * 8d0,      10.d0/
!         fractional energy of gamma Eg/Ek
       data kv/
     * 0.0, 0.050d0, 0.1d0,     0.15d0, 0.2d0, 0.25d0,
     * 0.3d0, 0.35d0, 0.4d0, 0.45d0, 0.5d0, 0.55d0,
     * 0.6d0, 0.65d0, 0.7d0, 0.75d0, 0.8d0, 0.85d0,
     * 0.9d0, 0.925d0, 0.95d0, 0.97d0, 0.99d0, 0.995d0,
     * 0.9990d0, 0.99950d0, 0.99990d0, 0.99995d0,
     * 0.99999d0, 1.0d0/

       if(first) then
          do i = 1, em
             Eklog(i) = log(Ek(i))
          enddo
          first = .false.
       endif

!       if(Ekin .lt. .001 .or.  intz .gt. 30) then
!          points = 3
!       elseif(Ekin .lt. .01 .or. intz .gt.20)  then
!          points = 4
!       elseif(intz .gt. 10) then
!          points = 5
!       else
!          points = 7
!       endif


       call kpolintp2(kv, 1, 0.d0, Eklog, 1, 0.d0,
     * dsdks, kmax, kmax, emax, points, 2, k, log(Ekin), ans, error)
       end function   epBrSeltzer



      subroutine epSeltzXsec(Z, dsdxs, kmax, emax)
      implicit none
!       get cross-section for atomic number Z

      integer Z   ! input 1~100

      integer kmax ! input. nsee below.  30
      integer emax ! input. see below.  <= 57
!        parameter (emax = 57, kmax = 30)
      real*8    dsdxs(kmax, emax)  ! ouput. read cross-section
                                   ! in (beta/Z)**2 v ds/dv in mb
      !   1~ kmax; index for Eg/(Ee-me) 
      !   1~ emax: index for Ee-me

      integer io
      character*4 name
      integer j, i

      data io/11/

      
      if(Z .lt. 1 .or. Z .gt. 100) then
         call cerrorMsg('Z is invalid for epSeltzXsec',0)
      endif
      
      if(Z .le. 9) then
         write(name, '("Z",i1)') Z
      elseif(Z .le. 99) then
         write(name, '("Z",i2)') Z
      else
         write(name, '("Z",i3)') Z
      endif

      open(io, file='Seltzer/'//name)

      do i = 1, emax
         do j = 1, kmax
            read(io, *)  dsdxs(j, i)
         enddo
         read(io, *)
      enddo
      close(io)
      end subroutine epSeltzXsec
