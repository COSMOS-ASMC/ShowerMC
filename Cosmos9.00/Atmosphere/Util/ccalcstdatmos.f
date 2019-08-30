!      Calculate standard atmospher constants.
!   stdatmos1.d is used to compute the consts.
!
!    This program outputs the following quantity, which may be saved in
!    stdatmos2.d which in turn is used by cstdatmos0.f to manage
!    atmosphere in the vertical direction.
!
!     height, temperature, pressure, rho, alfa, d0, dsum, H
!
!     all in mks unit.
!       dsum:  amount of atmospher above the node.
!       d0:    see the formula below.
!       alfa:  a in the formula below.
!
!     The scale height is approximated by a
!   number of stright lines as a function of height. The data in
!   stdatmos1.d gives height, temperatur, etc at  each nodal point.
!
!         The scale height, H, is expressed by H = H0 + a(z-z0)
!                                                = kT/mg      
!   in each region.
!   We neglect height dependence of gravitational accelleration g,
!   and the average mass of  air molecules, m.   
!   Since the data table gives T(z)= T0 + b(z-z0) at the nodal points, 
!   we can first get b,
!   and then a by a = dH/dz = k/mg * b.  H at a nodal point, z,  
!   is obtained as H(z) =kT(z)/mg.
!   The density is given by
!              rho = rho0 * (1+ a(z-z0)/H(z0))**(-1-1/a)      (a != 0)
!                  = rho0 * exp(- (z-z0)/H)            (a =0; hence H is const)
!
!   The gramage between  given heights, z1 and z2  is by
!
!              d = d0 *(fd(z1) - fd(z2))  where
!         
!            fd(z) = (1+ a(z-z0)/H(z0))**(-1/a)                   (a != 0)
!
!                  =  exp(-(z-z0)/H )                             (a = 0)
!
!   where  d0 = rho0*H(z0)
!
!    
      implicit none
#include  "Zmanagerp.h"
#include  "Zglobalc.h"
#include  "Zearth.h"

      integer nodes    ! real number of nodal points given in the table.
      integer icon, i,ios
      
      real*8 z(maxnodes), t(maxnodes), p(maxnodes), rho(maxnodes)
      real*8 scaleh(maxnodes), alfa(maxnodes), beta(maxnodes),
     *       d0(maxnodes)
      real*8 m, g, k, nuc, c2, kbymg, dtdz, dsum(maxnodes)
!      m: avrage molecule mass number of air
!      g: gravitational acceleration.
!      k: Boltzman's const.
!     c2: c**2
      parameter( m=14.5 *2, nuc =938.3e6, g=9.80665, c2=c*c,
     *          k= 8.617e-5, kbymg = k*c2/(m*nuc*g))

!        read basic data
      call copenf(TempDev, "stdatmos1%d", icon)
      if(icon .ne. 0) stop 9999
      call cskipComment(TempDev, icon)
      if(icon .ne. 0) stop 9999
!
      nodes = 0
      do while(1 .gt. 0 )
         read(TempDev, *, iostat=ios)
     *    z(nodes+1), t(nodes+1), p(nodes+1), rho(nodes+1)
         if(ios .ne. 0) goto 10
         nodes = nodes + 1
      enddo
 10   continue
      close(TempDev)
      write(*, *) " # of nodal points =", nodes, " kbymg=", kbymg

      do i =1, nodes-1
         dtdz = (t(i+1) - t(i))/ (z(i+1) - z(i))
         scaleh(i) = kbymg * t(i)                 ! at z(i)
         beta(i) = dtdz
         alfa(i) = kbymg * beta(i)
         d0(i) = rho(i) * scaleh(i)
!        write(*, *) z(i), d0(i), alfa(i), beta(i), scaleh(i)
      enddo

      dsum(nodes-1) = d0(nodes-1)

      do i = nodes-2, 1, -1
            if(alfa(i) .eq. 0.) then
!                exponential atmosphere
               dsum(i) = dsum(i+1) +
     *           d0(i)*(1.0- exp(-(z(i+1) - z(i))/scaleh(i)))
            else
               dsum(i) = dsum(i+1) +
     *            d0(i)*
     *           (1.0 -
     *             (1.0+alfa(i)*(z(i+1) -z(i))/scaleh(i))
     *              **(-1./alfa(i))
     *            )
            endif
      enddo
!          height, temp, pres, rho, alfa, d0, dsum, H
!       dsum is the grammage above the node. all in mks unit.
      write(*,*)"# height, temp, press, rho, alfa, d0, dsum, H"
      write(*,*) "#-----------------------------------------"
      do i = 1, nodes-1
         write(*,
     *   '(8g14.5)')
     *   sngl(z(i)),sngl(t(i)),sngl( p(i)), sngl(rho(i)),
     *   sngl(alfa(i)), sngl( d0(i)), sngl(dsum(i)), sngl(scaleh(i))
      enddo
!      do zx = 0., 11.d3, 1.d3
!         write(*, *) 
!     *    rho(1)*(1 + alfa(1)*(zx - 0.)/scaleh(1))**(-1.-1./alfa(1))
!      enddo
      end
