      program atmncOBPTCL
      implicit doubleprecision(a-h,o-z)
#include "Zcoord.h"
#include "Zcode.h"
      
      include '../src3/include/atmnc-particle-code.inc'
      logical,save::onlyfor1ry=.true.
      logical,save::onlyhadron=.true.
      integer:: lread(9)
      integer:: k1, k100
      real(8):: r1(5)  ! rx, ry, rz, t, r: r1(1:3)=same as comsmos
      real(8):: p1(5)  ! px,py,pz, E, m 
      real(8):: w(3)   ! direction cos.
      integer:: code, subc, charge
      integer:: sel
      type (coord):: a, b
      real(8):: p, r, cosz
      character(len=10) input

      if( iargc() == 1 ) then
         call getarg(1,input)
         read(input,*) sel
      else
         write(0,*)
     *   'Basic usage: ./translate$ARCH arg < Atmnc3-output '
         write(0,*) 'arg:'
         write(0,*)
     *   '     1==>selecting only hadronic particls suited for Cosmos'
         write(0,*)
     *   '         1ry. The content is sufficient for cmkInc2'
         write(0,*)
     *   '     2==>same as 1 but include e and g, if any'
         write(0,*)
     *   '     3==>show all particles; content is more than above'
         stop
      endif

      if(sel == 1 ) then
         ! default 
      elseif(sel == 2) then
         onlyhadron = .false.
      elseif(sel == 3) then
         onlyfor1ry = .false.
         onlyhadron = .false.
      else
         write(0,*) ' arg=', sel, ' invalid'
         stop
      endif   

      if(onlyfor1ry) then
         write(0,*) 
     *    'code subc charge   Wxyz  Et  Rxyz  llh  cosz'
         write(0,*) 
     *    '  1    2     3     4-6   7   8-10  11-13  14'
      else
         write(0,*) 
     *    'code subc charge k100  Pxyz  KE  Rxyz t  R   llh  cosz'
         write(0,*) 
     *    '  1    2     3    4    5 6 7  8  9-11 12 13  14-16  17'
      endif

      do while(.true.)  
         read(*,'(I4,5I9,3I11)',end=9999, err=1999) ! new file 
     &        (lread(i),i=1,9)
         call translate_aline(lread, k1, r1, p1, k100)
         call catmncTcos(k1, code, subc, charge)
         a%sys = 'xyz'
         a%r(1:3) = r1(1:3)
         call ctransCoord2('llh', a, b)
         p = sqrt( dot_product(p1(1:3), p1(1:3)))
         r = r1(5) !   =sqrt( dot_product(r1(1:3), r1(1:3)))
         cosz = dot_product(p1(1:3), r1(1:3))/p/r

         if(onlyfor1ry) then
            if( k100 < 3 .and. cosz < 0. .and.
     *         (code /= kneue .or. code /= kneumu) ) then
               if( (onlyhadron .and. code > kelec) .or.
     *              .not. onlyhadron )  then
                  w(1:3)= p1(1:3)/p
                  write(*, '(i3, i4, i3, 1p, 11g14.6)')
     *           code, subc, charge,  w(1:3), p1(4), r1(1:3),
     *           b%r(1:3), cosz
               endif
            endif
         else
            write(*,'(i3, i4, i3, i5, 1p, 13g14.6)')
     *           code, subc, charge,  k100, p1(1:3), p1(4)-p1(5),
     *           r1(1:5), b%r(1:3), cosz
         endif
         if(k100.ge.3) then ! primary, record here for all
!            call app_1ryratio(ik, ek, pickup)
!            crsamp = crsamp + 1
!            crsim  = crsim  + 1
         else
            if((k1 >= knue) .and. (k1 <= knumubar)) then
               ! nutrino
            else
               ! r1(5) is the crossing radius
            end if
         end if                 ! neutrino or primary
 1999    continue               ! anyway move to next file
      enddo

 9999 continue
      end      program atmncOBPTCL

      subroutine translate_aline(lread, k1, r1, p1, k100)
      implicit doubleprecision (a-h, o-z)
! masses of particles
      include '../src3/include/atmnc-particle-mass2.inc'

      dimension r1(5), p1(5)
      dimension lread(9)

      kk   = lread(1)
      k100 = int(kk/100)
      k1   = kk - 100 * k100
!      write(*,*) lread(1), k100, k1

      pp = 10.d0**(1.d-6*lread(2))
      rr = 0.
      do k     = 1, 3
         r1(k) =        lread(2+k)
         p1(k) = 1.d-9 * lread(6+k)*pp
         rr    = rr + r1(k)**2
      enddo
      rr = sqrt(rr)

      r1(4) = 10.d0**(1.d-6*real(lread(6)))
      r1(5) = rr

      p1(5) = am(k1)
      p1(4) = sqrt(pp**2 + p1(5)**2)
      end      subroutine translate_aline
