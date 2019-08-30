!c      test cmucap
!      real*8 capr
!       call cmucap(16,  8,    capr)
!       write(*,*) ' capr=',capr
!       end
!  compute muon capture rate /sec
!    mu- + p --> n + neu(mu)
!    mu- + (z,n)--->(z+1, n-1) + neu(mu)
!   at  z=9, decay rate and capture rate compeat
!
         subroutine cmucap(a, z, capr)
!          a: integer. input.  mass number of the medium
!          z: integer. input.  atomic number of the medium
!       capr: real*8. output. capture rate/sec
!
         implicit none
         integer a, z
         real*8 capr
!
!        integer nn
!        parameter (nn=16)
!c       real mass(nn)/
!    *   9.00, 12.01, 16.00, 24.32, 28.11, 32.08, 40.10,
!    *  51.00, 58.77,
!    *   63.62, 95.98, 112.52, 184.0, 207.8, 209.0,
!    *  238.0/
!        real rho(nn)/
!    *   4.528e-6, 2.026e-5, 3.786e-5, 2.478e-4, 4.265e-4,
!    *   6.612e-4, 1.290e-3, 2.016e-3, 3.462e-3,
!    *    3.789e-3,9.189e-3, 1.204e-2, 2.190e-2,
!    *    2.592e-2,
!    *   2.542e-2, 2.832e-2/
!
!     reference:1)h. primakoff. rev. mod. phys. 31(1959)802
!               2)j.c.sens, r.a.swanson, v.l.telegdi, d.d.yovanovitch
!                 p.r. 107(1957)1464.
!               3)j.c.sens, p.r. 113 (1959)679 (refer this one)
!          capture rate for o:  125000/sec
!                           n:   86000/sec
!                          ar: 1530000/sec
!                      s.rock: 360000/sec
         real*8 x, rho
!
         if(a .gt. 8) then
!              for the case of air, treat specially.
            if(a .eq. 14 .and. z .eq. 7) then
                 capr=0.86d5
            elseif(a .eq. 16 .and.  z .eq. 8) then
!                   experimental value of o is not reliable.
!                   use (zeff(o)/zeef(c))**4 * capr(c)
!                   as (7.47/5.75)**4 * 0.44e5=125000.
                 capr=125000.d0
            else
                 x=log10(float(a))
!                   nuclear density (protons/fermi**3)
                 rho=10.0**( (-1.37*x +7.23)*x -11.0)
                 capr=rho *
     *           (- 21.8/0.07 *( (a-z)/2/float(a)-0.24) + 24.05)*1.e8
            endif
         else
             capr=169. * z**4 * (1. - 3.15 * (a-z)/2/float(a) )
         endif
      end
