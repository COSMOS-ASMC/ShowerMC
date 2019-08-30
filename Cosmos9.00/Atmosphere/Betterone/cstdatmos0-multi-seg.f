!    Programs contained in this file  manage the atmosphere in 
!   the vertical direction.
!   stdatmos2.d is used to describe the atmosphere. > 11km
!    
!      cstdatmos0:  to read data for the atmospeher
!      cvh2den   :  vertical height to denstiy
!      cvh2thick :  vertical height to thickness of air
!      cvthick2h :  vertical thicknes to height.
!      cvh2denp  :  d rho/ dz
!      cvh2den2p  :  d (d rho/ dz)/dz
!

!    They follow the folloing formulas.
!
!          rho = rho0 * (1+ a(z-z0)/H(z0))**(-1-1/a)          (a != 0)
!              = rho0 * exp(- (z-z0)/H)           (a =0; hence H is const)

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
!   Other formula is used  at < 11 km. I.e.,
!
!         rho = rho00*( (ha - z)/hl )** pw
!
!
      

!    
      subroutine cstdatmos0
      implicit none
#include  "Zmanager.h"
#include  "Zearth.h"
#include  "Zstdatmos.h"

      integer icon,  ios, i
      external blkstdatmos
c----      include 'Zstdatmosf.h'     !  rho=rho(h) function
#include  "Zstdatmosf.h"

      if(first) then
         call copenf(TempDev, "stdatmos2%d", icon)
         if(icon .ne. 0) stop 9999
         call cskipComment(TempDev, icon)
         if(icon .ne. 0) stop 9999

         nodes = 0
         do while( .true. )
            read(TempDev, *, iostat=ios)
     *           znode(nodes+1), tnode(nodes+1), pnode(nodes+1), 
     *           rho(nodes+1),  alfa(nodes+1), d0(nodes+1),
     *           dsum(nodes+1),  scaleh(nodes+1)
            if(ios .ne. 0) goto 10
            if(nodes .ge. maxnodes) then
               write(*,*) 'numbr of nodes for atmosphere > ', maxnodes
               stop 9999
            endif
            nodes = nodes + 1
!           write(0,*) znode(nodes), alfa(nodes)
         enddo
 10      continue
!        write(*, *) " # of nodal points =", nodes

         do i=2, nodes
            if(alfa(i-1) .ne. 0.) then
               fd1i(i) = fd1(znode(i), alfa(i-1), znode(i-1), 
     *                   scaleh(i-1))
               rhop(i) =  rho(i-1) *(-1.d0 -1.d0/alfa(i-1)) *
     *              alfa(i-1)/scaleh(i-1) 
               pwp(i) =-2.d0-1.d0/alfa(i-1)
            else
               fd0i(i) = fd0(znode(i), znode(i-1), scaleh(i-1) )
            endif
         enddo
         fd3 = ((ha-znode(2))/hl)**(pw+1.d0)
         first = .false.
      endif
      end
!     *************************    to be updated 
      block data blkstdatmos
!     *************************
       implicit none
c----       include 'Zearth.h'
#include  "Zearth.h"
c----       include 'Zstdatmos.h'
#include  "Zstdatmos.h"
       data first/.true./
       data ha/45.3425d3/, hl/11.865d3/,
     *     pw /4.36815d0/, rho00/3.50618d-3/
      end
!     *********************************   
      real*8 function  cvh2den(z)
!     *********************************   
!         gives density of air at a given vertical height.
!      h: vertical height in m
!      function value:  density in kg/m3.
!
      implicit none
c----      include 'Zearth.h'
#include  "Zearth.h"
c----      include 'Zstdatmos.h'
#include  "Zstdatmos.h"
!
      integer i

      real*8 z

      if(first) call cstdatmos0

!          check if z < 11.2km
      if(z .lt. znode(2)) then
         cvh2den = rho00 * ( (ha - z)/hl )**pw
      else
         do i = 3, nodes
            if(z .lt. znode(i) .or. i .eq. nodes ) then
               if(alfa(i-1) .ne. 0.)then
                  if(z .lt. znode(i-1) ) then
                     write(0,*) ' z=',z, ' znode =', znode(i-1)
                  endif
                  cvh2den = 
     *             rho(i-1) * 
     *             (1.d0 + alfa(i-1)* 
     *             (z-znode(i-1))/scaleh(i-1) )**(-1.0-1.d0/alfa(i-1)) 
                   goto 10
               else
                  cvh2den =
     *            rho(i-1) * exp(- (z-znode(i-1))/scaleh(i-1))         
                  goto 10
               endif
            endif
         enddo
 10      continue
      endif
      end
!     *************************
      real*8 function cvh2thick(z)
!         vertical height to thickness of air above that point.
!      z:  input. vertical height in m
!  function value. output.  thickness in kg/m2.
!
!
      implicit none
c----      include 'Zearth.h'
#include  "Zearth.h"
c----      include 'Zstdatmos.h'
#include  "Zstdatmos.h"
!
      integer i

      real*8 z

c----      include 'Zstdatmosf.h'
#include  "Zstdatmosf.h"


      if(first) call cstdatmos0

!          check if z < 11.2km
      if(z .lt. znode(2)) then
           cvh2thick =hl* rho00 /(pw+1.d0)* (
     *         ( (ha-z)/hl)**(pw+1.d0) -
     *         fd3 )
     *      +  dsum(2)
!              where fd3 is
!     *           ((ha-znode(2))/hl)**(pw+1.d0) )
      else
!          The gramage between  given heights, z1 and z2  is by
!
!              d = d0 *(fd(z1) - fd(z2))  where
!         
!            fd(z) = (1+ a(z-z0)/H(z0))**(-1/a)                   (a != 0)
!
!                  =  exp(-(z-z0)/H )                             (a = 0)
         do i = 3, nodes-1
            if(z .lt. znode(i) .or. i .eq. nodes-1 ) then
               if(alfa(i-1) .ne. 0.)then
                   cvh2thick =  d0(i-1) * (
     *              fd1( z, alfa(i-1), znode(i-1), scaleh(i-1) )
     *              - fd1i(i)
     *             )   + dsum(i) 

!                      where fd1i(i) is
!     *              fd1(znode(i), alfa(i-1),
!     *              znode(i-1), scaleh(i-1) )

                   goto 10
               else
                  cvh2thick = d0(i-1) * (
     *             fd0(z, znode(i-1), scaleh(i-1) )
     *             - fd0i(i)               
     *             ) + dsum(i)
!                      where fd0i is
!     *            fd0(znode(i), znode(i-1), scaleh(i-1) )
!     *            )
                  goto 10
               endif
            endif
         enddo
 10      continue
      endif
      end
!     *************************
      real*8 function cvthick2h(t)
!         vertical thickness to vertical height conversion
!     t:  input. vertical tikness in kg/m2
!  function value. output.  height in m.
!
!
      implicit none
c----      include 'Zearth.h'
#include  "Zearth.h"
c----      include 'Zstdatmos.h'
#include  "Zstdatmos.h"
!
      integer i

      real*8  t,  temp


      if(first) call cstdatmos0

!          check if z < 11.2km
      if(t .gt. dsum(2)) then
!           t =hl* rho00 /(pw+1.d0)* (
!     *         ( (ha-z)/hl)**(pw+1.d0) -
!     *           ((ha-znode(2))/hl)**(pw+1.d0) )
!     *      +  dsum(2)
!           solve above eq.
            temp =( t - dsum(2))/hl/rho00 * (pw+1.d0)
     *          +  fd3
            cvthick2h =ha -  hl* temp**(1.d0/(pw+1.d0))
      else
!          The gramage between  given heights, z1 and z2  is by
!
!              d = d0 *(fd(z1) - fd(z2))  where
!         
!            fd(z) = (1+ a(z-z0)/H(z0))**(-1/a)                   (a != 0)
!
!                  =  exp(-(z-z0)/H )                             (a = 0)
         do i = 3, nodes-1
            if(t .gt. dsum(i) .or. i .eq. nodes-1 ) then
               if(alfa(i-1) .ne. 0.)then
!                   cvh2thick =  d0(i-1) * (
!     *              fd1( z, alfa(i-1), znode(i-1), scaleh(i-1) )
!     *             -  fd1(znode(i), alfa(i-1),
!     *              znode(i-1), scaleh(i-1) )
!     *             )   + dsum(i) 
!                    solve above eq.
                   temp =( t -dsum(i))/d0(i-1)
     *              +  fd1i(i)
                   cvthick2h = (temp**(-alfa(i-1)) -1.d0)*
     *               scaleh(i-1)/alfa(i-1) + znode(i-1)
                   goto 10
               else
!                  cvh2thick = d0(i-1) * (
!     *             fd0(z, znode(i-1), scaleh(i-1) )
!     *            -fd0(znode(i), znode(i-1), scaleh(i-1) )
!     *            )  + dsum(i)
!                    solve above eq.
                   temp = (t -dsum(i)) /d0(i-1) +
     *             + fd0i(i)
!                   temp =  exp(-(z-z0)/H )    
                   cvthick2h =znode(i-1)- log(temp)*scaleh(i-1)
                  goto 10
               endif
            endif
         enddo
 10      continue
      endif
      end
!     *********************************   
      real*8 function  cvh2denp(z)
!     *********************************   
!         gives derivative of the density of
!     air at a given vertical height.
!      h: vertical height in m
!      function value: drho/dz density in kg/m4
!
      implicit none
c----      include 'Zearth.h'
#include  "Zearth.h"
c----      include 'Zstdatmos.h'
#include  "Zstdatmos.h"
!
      integer i

      real*8 z

      if(first) call cstdatmos0

!          check if z < 11.2km
      if(z .lt. znode(2)) then
         cvh2denp =- rho00 * pw/hl* ( (ha - z)/hl )**(pw-1.d0)
      else
         do i = 3, nodes
            if(z .lt. znode(i) .or. i .eq. nodes ) then
               if(alfa(i-1) .ne. 0.)then
                  cvh2denp = 
!       *             rho(i-1) *(-1.d0 -1.d0/alfa(i-1)) *
!       *             alfa(i-1)/scaleh(i-1) *     = rhop(i)
     *               rhop(i) *
     *             (1.d0 + alfa(i-1)* 
!     *             (z-znode(i-1))/scaleh(i-1) )**(-2.d0-1.d0/alfa(i-1)) 
     *             (z-znode(i-1))/scaleh(i-1) )**pwp(i)
                   goto 10
               else
                  cvh2denp =
     *            -rho(i-1) * exp(- (z-znode(i-1))/scaleh(i-1))         
     *            /scaleh(i-1)
                  goto 10
               endif
            endif
         enddo
 10      continue
      endif
      end
!     *********************************   
      real*8 function  cvh2den2p(z)
!     *********************************   
!         gives double derivative of the density of
!     air at a given vertical height.
!      h: vertical height in m
!      function value: d (drho/dz)/dz density in kg/m5
!
      implicit none
c----      include 'Zearth.h'
#include  "Zearth.h"
c----      include 'Zstdatmos.h'
#include  "Zstdatmos.h"
!
      integer i

      real*8 z

      if(first) call cstdatmos0

!          check if z < 11.2km
      if(z .lt. znode(2)) then
         cvh2den2p = rho00 * pw*(pw-1.d0)/hl/hl
     *          * ( (ha - z)/hl )**(pw-2.d0)
      else
         do i = 3, nodes
            if(z .lt. znode(i) .or. i .eq. nodes ) then
               if(alfa(i-1) .ne. 0.)then
                  cvh2den2p = 
     *               rhop(i) * pwp(i) *alfa(i-1)/scaleh(i-1)*
     *             (1.d0 + alfa(i-1)* 
     *             (z-znode(i-1))/scaleh(i-1) )**(pwp(i)-1.d0)
                   goto 10
               else
                  cvh2den2p =
     *            rho(i-1) * exp(- (z-znode(i-1))/scaleh(i-1))         
     *            /scaleh(i-1)**2
                  goto 10
               endif
            endif
         enddo
 10      continue
      endif
      end

