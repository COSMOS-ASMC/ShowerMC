      subroutine ciniSegAtmos
!        read the atmosphere data
!
      implicit none
#include  "Zmanager.h"
#include  "Zmanagerp.h"
#include  "Zearth.h"
#include  "Zatmos.h"

      integer icon,  ios, i
      character*100 msg

      call copenf(TempDev, AtmosFile, icon)
      if(icon .ne. 0) stop 9999
      call cskipComment(TempDev, icon)
      if(icon .ne. 0) stop 9999

      nodes = 0
      mostz = 0
      do while( .true. )
         read(TempDev, *, iostat=ios)
     *      atmos%z(nodes+1), atmos%T(nodes+1), atmos%P(nodes+1), 
     *      atmos%rho(nodes+1),  atmos%a(nodes+1), atmos%d0(nodes+1),
     *      atmos%cumd(nodes+1),  atmos%H(nodes+1)
         if(ios .ne. 0) goto 10
         if(nodes .ge. maxnodes) then
            write(msg,*) 'numbr of nodes for atmosphere > ', maxnodes
            call cerrorMsg(msg, 0)
         endif
         nodes = nodes + 1
         if(mostz .eq. 0 .and.   Znode(nodes) .gt. 2900.) then
            mostz = nodes
         endif
!        write(0,*) Znode(nodes), Anode(nodes)
      enddo
 10   continue
!     write(*, *) " # of nodal points =", nodes, ' mostz', mostz

      do i=1, nodes-1
         if(Anode(i) .ne. 0.) then
            fd1i(i) = fd1(Znode(i+1), Anode(i), Znode(i), 
     *           Scalehnode(i))
            Rhopnode(i) =  Rhonode(i) *(-1.d0 -1.d0/Anode(i)) *
     *           Anode(i)/Scalehnode(i) 
            pwp(i) =-2.d0-1.d0/Anode(i)
         else
            fd0i(i) = fd0( Znode(i+1), Znode(i), Scalehnode(i) )
         endif
      enddo
      endif
      end
!     *********************************   
      real*8 function  cvh2den(z)
!     *********************************   
!         gives density of air at a given vertical height.

!
      implicit none
#include  "Zatmos.h"

      real*8 z     ! input.  vertical height in m
!                   function value:  density in kg/m3.
      integer i 

      if( z .ge. Znode(mostz) .and. z .lt. Znode(mostz+1)) then
         i = mostz
      elseif( z .ge Znode(mostz + 1) .and. z .lt. Znode(mostz+2) then
         i = mostz + 1
      else
         do i = mostz, 2, -1
            if(z .ge. Znode(i-1) .and. z .lt. Znode(i)) then
               goto 100
            endif
         enddo
         call kdwhereis(z, nodes, Znode, 1, i)
         if(i .eq. 0) then
            call cerrorMsg('height becomes too small ', 0)
         endif
      endif
 100  continue
      if(Anode(i) .ne. 0.)then
         cvh2den = 
     *        Rhonode(i) * 
     *       (1.0 + Anode(i)* 
     *          (z-Znode(i))/Scalehnode(i) )**(-1.0-1.d0/Anode(i)) 
      else
         cvh2den =
     *        Rhonode(i) * exp(- (z-Znode(i))/Scalehnode(i))         
      endif
      end
!     *************************
      real*8 function cvh2thick(z)
!           vertical height to thickness of air above that point.
!      z:  input. vertical height in m
!  function value. output.  thickness in kg/m2.
!
!
      implicit none
#include  "Zearth.h"
#include  "Zstdatmos.h"
#include  "Zatmos.h"

      real*8 z        ! input
      integer i

#include  "Zstdatmosf.h"


!          The gramage between  given heights, z1 and z2  is by
!
!              d = D0node *(fd(z1) - fd(z2))  where
!         
!            fd(z) = (1+ a(z-z0)/H(z0))**(-1/a)                   (a != 0)
!
!                  =  exp(-(z-z0)/H )                             (a = 0)

      if( z .ge. Znode(mostz) .and. z .lt. Znode(mostz+1)) then
         i = mostz
      elseif( z .ge Znode(mostz + 1) .and. z .lt. Znode(mostz+2) then
         i = mostz + 1
      else
         do i = mostz, 2, -1
            if(z .ge. Znode(i-1) .and. z .lt. Znode(i)) then
               goto 100
            endif
         enddo
         call kdwhereis(z, nodes, Znode, 1, i)
         if(i .eq. 0) then
            call cerrorMsg('height becomes too small ', 0)
         endif
      endif
 100  continue

      if(Anode(i) .ne. 0.)then
         cvh2thick = Dsumnode(i) -
     *     D0node(i) * (
     *     1.0- fd1(z, Anode(i), Znode(i), Scalehnode(i)) 
     *                 )
      else
         cvh2thick = Dsumnode(i) +
     *     D0node(i) * (1. -
     *       fD0(z, Znode(i), Scalehnode(i) )
     *                 )
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
#include  "Zatmos.h"

      integer i

      real*8  t,  temp


      if( z .ge. Znode(mostz) .and. z .lt. Znode(mostz+1)) then
         i = mostz
      elseif( z .ge Znode(mostz + 1) .and. z .lt. Znode(mostz+2) then
         i = mostz + 1
      else
         do i = mostz, 2, -1
            if(z .ge. Znode(i-1) .and. z .lt. Znode(i)) then
               goto 100
            endif
         enddo
         call kdwhereis(z, nodes, Znode, 1, i)
         if(i .eq. 0) then
            call cerrorMsg('height becomes too small ', 0)
         endif
      endif
 100  continue

      if(Anode(i) .ne. 0. )then
!         cvh2thick = Dsumnode(i) -
!     *     D0node(i) * (
!     *     1.0- fd1(z, Anode(i), Znode(i), Scalehnode(i)) 
!     *                 )
!                    solve above eq.
         temp =1.0 -  (Dsumnode(i) - t )/D0node(i)

         cvthick2h = (1.0 -  temp**(-Anode(i)))*
     *               Scalehnode(i)/Anode(i) + Znode(i)
      else
!          t = Dsumnode(i) +
!     *     D0node(i) * (1. -
!     *       fD0(z, Znode(i), Scalehnode(i) )
!     *                 )
!            solve above eq.
         temp =1.0- (t -Dsumnode(i)) /D0node(i) 
!                   temp =  exp(-(z-z0)/H )    
         cvthick2h =Znode(i)- log(temp)*Scalehnode(i)
      endif
      end
!     *********************************   
      real*8 function  cvh2denp(z)
!     *********************************   
!         gives derivative of the density of
!     air at a given vertical height.
!      h: vertical height in m
!      function value: dRhonode/dz density in kg/m4
!
      implicit none
c----      include 'Zearth.h'
#include  "Zearth.h"
c----      include 'Zstdatmos.h'
#include  "Zstdatmos.h"
!
      integer i

      real*8 z


      if( z .ge. Znode(mostz) .and. z .lt. Znode(mostz+1)) then
         i = mostz
      elseif( z .ge Znode(mostz + 1) .and. z .lt. Znode(mostz+2) then
         i = mostz + 1
      else
         do i = mostz, 2, -1
            if(z .ge. Znode(i-1) .and. z .lt. Znode(i)) then
               goto 100
            endif
         enddo
         call kdwhereis(z, nodes, Znode, 1, i)
         if(i .eq. 0) then
            call cerrorMsg('height becomes too small ', 0)
         endif
      endif
 100  continue

      if(Anode(i) .ne. 0.)then
         cvh2denp = 
!       *             Rhonode(i-1) *(-1.D0node -1.D0node/Anode(i-1)) *
!       *             Anode(i-1)/Scalehnode(i-1) *     = Rhonodep(i)
     *        Rhopnode(i) *
     *        (1.d0 + Anode(i)* 
!     *             (z-Znode(i))/Scalehnode(i) )**(-2.d0-1.d0/Anode(i)) 
     *       (z-Znode(i))/Scalehnode(i) )**pwp(i)
      else
         cvh2denp =
     *        -Rhonode(i) * exp(- (z-Znode(i))/Scalehnode(i))         
     *        /Scalehnode(i)
      endif
      end
!     *********************************   
      real*8 function  cvh2den2p(z)
!     *********************************   
!         gives double derivative of the density of
!     air at a given vertical height.
!      h: vertical height in m
!      function value: d (dRhonode/dz)/dz density in kg/m5
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
      if(z .lt. Znode(2)) then
         cvh2den2p = Rhonode00 * pw*(pw-1.D0node)/hl/hl
     *          * ( (ha - z)/hl )**(pw-2.D0node)
      else
         do i = 3, nodes
            if(z .lt. Znode(i) .or. i .eq. nodes ) then
               if(Anode(i-1) .ne. 0.)then
                  cvh2den2p = 
     *               Rhonodep(i) * pwp(i) *Anode(i-1)/Scalehnode(i-1)*
     *             (1.D0node + Anode(i-1)* 
     *             (z-Znode(i-1))/Scalehnode(i-1) )**(pwp(i)-1.D0node)
                   goto 10
               else
                  cvh2den2p =
     *            Rhonode(i-1) * exp(- (z-Znode(i-1))/Scalehnode(i-1))         
     *            /Scalehnode(i-1)**2
                  goto 10
               endif
            endif
         enddo
 10      continue
      endif
      end

