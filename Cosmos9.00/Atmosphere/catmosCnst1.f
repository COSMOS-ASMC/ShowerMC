      subroutine catmosCnst1

!   catmosCnst1:  compute some basic constants
!      a,  b,  d0, cum, H
! where
!     a, b; const shown below
!       d0: see the formula below.
!     cumd: amount of atmospher above the node.
!        H: scale height at the node
!     all in mks unit.
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
!
!   The density off a nodal point is given by
!
!              rho = rho0 * (1+ a(z-z0)/H(z0))**(-1-1/a)      (a != 0)
!                  = rho0 * exp(- (z-z0)/H)            (a =0; hence H is const)
!
!   (We employ H(z0) as the scale height in the segment)
!   The amount of air between  given heights, z1 and z2  is by
!
!              d = d0 *(fd(z1) - fd(z2))  where
!         
!            fd(z) = (1+ a(z-z0)/H(z0))**(-1/a)                   (a != 0)
!
!                  =  exp(-(z-z0)/H )                             (a = 0)
!
!   where  d0 = rho0*H(z0)
!
!   If z1=z0,  d becomes
!
!              d= d0 ( 1 - fd(z2))
!
!

!
      use modAtmosDef  
      implicit none

#include  "Zglobalc.h"
#include  "ZmediaLoft.h"

      integer i

      

      real*8 m,  k, nuc, c2, kbymg, dtdz
!      m: avrage molecule mass number of air
!     gaccel: = old g; gravitational acceleration.  m/sec2
!        should be defined in modAtmosDef.
!      k: Boltzman's const.
!     c2: c**2
!      parameter( m=14.5d0 *2, nuc =938.3d6, g=9.80665d0, c2=c*c,
!     *          k= 8.617d-5, kbymg = k*c2/(m*nuc*g))
      parameter(  nuc = 938.3d6,  c2=c*c,
     *          k= 8.617d-5)



      integer nodes, mno

#include  "Zstdatmosf.h"

      real(8),external:: cvh2thick

      nodes = atmos%nodes
      do i = 1, nodes-1
         mno= atmos%node2mediaNo(i)
!          use only for gas; but compute for solid too to avoid
!           undeef when printing.
         m = Media(mno)%A * Media(mno)%sumNo
         kbymg = k*c2/(m*nuc*gaccel)
         dtdz =
     *      (atmos%T(i+1) - atmos%T(i))/ (atmos%z(i+1) - atmos%z(i))

         atmos%H(i) = kbymg * atmos%T(i) ! at atmos%z(i)
         atmos%b(i) = dtdz
         atmos%a(i) = kbymg * atmos%b(i)
         atmos%d0(i) = atmos%rho(i) * atmos%H(i)

         if( Media(mno)%gasF == 0 ) then
           !    non gas; 
           ! except rho, all above will not be used. 
            atmos%b(i) = 0.  ! but put some
            atmos%a(i) = 0.  !  here
            atmos%rho(i) = Media(mno)%rho * atmos%rhoc(i)*1.0d3 !g/cm3-->kg/m3
!            Meias(mno)%rhoc is kept 1.0. later it will be changed
!     for consistency with other cases
            atmos%rhoc(i) = -atmos%rhoc(i) ! non -gas flag.  gas will be 1.0
         endif

!///////////
!        write(0, *)
!     * atmos.z(i), atmos.d0(i), atmos.a(i), atmos.b(i), atmos.H(i)
!////////
      enddo

!            top node must be gas.
      if( atmos%rhoc(nodes) < 0. ) then
         write(0,*) ' the top node matter is =', atmos%matter(nodes)
         write(0,*) ' it should be gas but not so '
         write(0,*) ' For a Moon like object, put "sp" at the top'
         stop
      else
         atmos%cumd(nodes) =  atmos%rho(nodes) * Hinf ! put  very small amout
      endif

      do i = nodes-1, 1, -1
         if( atmos%rhoc(i) > 0. ) then         
            if(atmos%a(i) .eq. 0.) then
!                 exponential atmosphere
               atmos%cumd(i) = atmos%cumd(i+1) +
     *           atmos%d0(i)*
     *         (1.0- fd0(atmos%z(i+1), atmos%z(i), atmos%H(i)))
            else
               atmos%cumd(i) = atmos%cumd(i+1) +
     *              atmos%d0(i)*
     *              (1.0 -
     *         fd1(atmos%z(i+1), atmos%a(i),atmos%z(i), atmos%H(i) ) )
            endif
         else
       !    const rho non gas
            atmos%cumd(i) =  atmos%cumd(i+1)  +
     *             atmos%rho(i) * ( atmos%z(i+1) - atmos%z(i) )
         endif
      enddo

!        simply take log10      
      do i = 1, nodes
         atmos%logP(i) = log(atmos%P(i))
         atmos%logrho(i) = log(atmos%rho(i))
         atmos%logcumd(i) = log(atmos%cumd(i))
         atmos%logcumdi(atmos%nodes-i+1) = atmos%logcumd(i)
         atmos%logrhoi(atmos%nodes-i+1) = atmos%logrho(i)
         atmos%zi(atmos%nodes-i+1) = atmos%z(i)
      enddo

      AlmostVacT = cvh2thick( AlmostVacH )
      
      end

      subroutine catmosCnst2  ! would  not be used now
!     compute c-spline coef. for later use
      use modAtmosDef
      implicit none
!  #include "Zatmos.h"

      integer  nodes
      
      nodes = atmos%nodes
!                           height--> rho
      call kcsplCoef(atmos%z, atmos%logrho, nodes, atmos%coefh2r,
     *      maxnodes)
!                           height--> log(depth)
      call kcsplCoef(atmos%z, atmos%logcumd, nodes, atmos%coefh2d,
     *      maxnodes)
!                           h--> log(P)
!      call kcsplCoef(atmos.z, atmos.logP, nodes, atmos.coefh2P,
!     *      maxnodes)
!                           h--> scale H
      call kcsplCoef(atmos%z, atmos%H, nodes-1, atmos%coefh2H,
     *      maxnodes)
!                           h--> T
      call kcsplCoef(atmos%z, atmos%T, nodes, atmos%coefh2T,
     *      maxnodes)
!                           log(rho)--> h
!      call kcsplCoef(atmos.logrho, atmos.z, nodes, atmos.coefr2h,
!     *      maxnodes)
!                           log(depth)--> h
      call kcsplCoef(atmos%logcumdi, atmos%zi, nodes, atmos%coefd2h,
     *      maxnodes)

!               log(detph) --> log(rho)
      call  kcsplCoef(atmos%logcumdi, atmos%logrhoi, nodes,
     *      atmos%coefd2r, maxnodes)
!                           log(P) --> h 
!      call kcsplCoef(atmos.logP, atmos.z, nodes, atmos.coefP2h,
!     *      maxnodes)


      end
