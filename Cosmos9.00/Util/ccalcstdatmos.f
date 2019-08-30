!      Calculate standard atmospher constants.
!   stdatmos1.d is used to compute the consts.
!
!    This program outputs the following quantity, which may be saved in
!    stdatmos2.d which in turn is used by cstdatmos0.f to manage
!    atmosphere in the vertical direction.
!
!     height, temperature, pressure, rho, Anode, d0, dsum, H
!
!     all in mks unit.
!       dsum:  amount of atmospher above the node.
!       d0:    see the formula below.
!       Anode:  a in the formula below.
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
!   (We employ H(z0) as the scale height in the segment)
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
#include  "Zatmos.h"

      integer icon, i,ios

      

      real*8 m, g, k, nuc, c2, kbymg, dtdz
!      m: avrage molecule mass number of air
!      g: gravitational acceleration.
!      k: Boltzman's const.
!     c2: c**2
      parameter( m=14.5 *2, nuc =938.3e6, g=9.80665, c2=c*c,
     *          k= 8.617e-5, kbymg = k*c2/(m*nuc*g))

#include  "Zstdatmosf.h"


!        read basic data
      call copenf(TempDev, "../Data/Atmos/stdatmos1.d", icon)
      if(icon .ne. 0) stop 9999
      call cskipComment(TempDev, icon)
      if(icon .ne. 0) stop 9999
!
      nodes = 0
      do while(1 .gt. 0 )
         read(TempDev, *, iostat=ios)
     *   Znode(nodes+1), Tnode(nodes+1), Pnode(nodes+1),
     *   Rhonode(nodes+1)
         if(ios .ne. 0) goto 10
         nodes = nodes + 1
      enddo
 10   continue
      close(TempDev)
      write(*, *) "# of nodal points =", nodes, " kbymg=", kbymg

      do i = 1, nodes-1
         dtdz = (Tnode(i+1) - Tnode(i))/ (Znode(i+1) - Znode(i))
         Scalehnode(i) = kbymg * Tnode(i)                 ! at Znode(i)
         Bnode(i) = dtdz
         Anode(i) = kbymg * Bnode(i)
         D0node(i) = Rhonode(i) * Scalehnode(i)
!        write(*, *)
!    * Znode(i), D0node(i), Anode(i), Bnode(i), Scalehnode(i)
      enddo

      Dsumnode(nodes) =  0.

      do i = nodes-1, 1, -1
         if(Anode(i) .eq. 0.) then
!              exponential atmosphere
            Dsumnode(i) = Dsumnode(i+1) +
     *       D0node(i)*
!     *       (1.0- exp(-(Znode(i+1) - Znode(i))/Scalehnode(i)))
     *       (1.0- fd0(Znode(i+1), Znode(i), Scalehnode(i)))
         else
            Dsumnode(i) = Dsumnode(i+1) +
     *        D0node(i)*
     *        (1.0 -
     *          fd1(Znode(i+1), Anode(i),Znode(i), Scalehnode(i) ) )

!     *        (1.0+Anode(i)*(Znode(i+1) -Znode(i))/Scalehnode(i))
!     *        **(-1./Anode(i))
!     *         )
         endif
      enddo
!          height, temp, pres, Rhonode, Anode, D0node, Dsumnode, H
!       Dsumnode is the grammage above the node. all in mks unit.
      write(*,'(a)')
     * "# Height    Temp     Press   Rhonode',
     * '     Anode   D0node  Dsumnode    H"
      write(*,'(a)') "#-----------------------------------------"
      do i = 1, nodes-1
         write(*,
     *   '(8g14.5)')
     *   sngl(Znode(i)), sngl(Tnode(i)), sngl( Pnode(i)),
     *   sngl(Rhonode(i)), sngl(Anode(i)), sngl( D0node(i)),
     *   sngl(Dsumnode(i)), sngl(Scalehnode(i))
      enddo
      end


