!  This is to express atmosphere using linear interpolation.
!  cspline is not used because it generates some odd behaviour
!  
!
!      -------------------------------------
       real*8 function cvh2den(z)
!     --------------------------vertical height to density
       use modAtmosDef
       implicit none
! #include  "Zatmos.h"
       real*8 z  ! input. m
!       real*8 zsave
       real*8  ans
!       data zsave/-1.d30/
!       save zsave, ans
        integer i
        real*8 a

!       if(z .ne. zsave) then
          if( z .gt. atmos%z(atmos%nodes) ) then
             ans = atmos%rho(atmos%nodes)*
     *            exp(-(z-atmos%z(atmos%nodes))/Hinf)
          elseif(z .lt. atmos%z(1)) then
             ans = atmos%rho(1)*
     *          exp( (atmos%z(1)-z)/atmos%H(1) )
          else
             call kdwhereis(z, atmos%nodes, atmos%z, 1, i)
             a = atmos%a(i)
             if(a .ne. 0.d0) then
                  ans =  atmos%rho(i)*
     *            (1+ a*(z-atmos%z(i))/atmos%H(i))**(-1.0d0-1.d0/a)
	     else
                  ans =
     *             atmos%rho(i) * exp(- (z-atmos%z(i))/atmos%H(i))
            endif
          endif 
!          zsave = z
!       endif
       cvh2den = ans
       end
!     ----------------------------------
      real*8 function cvh2temp(z)
      use modAtmosDef
      implicit none
! #include "Zatmos.h"
!          vettical height to temperatur (Kelvin)

      real*8 z   ! input.  vertical height in m
!        output is temperature of the atmospher in Kelvin

      real*8 ans
      integer i

      if( z .gt. atmos%z(atmos%nodes) ) then
         ans = atmos%T(atmos%nodes)
      elseif(z .lt. atmos%z(1)) then
         ans = atmos%T(1) + atmos%b(1)*(z - atmos%z(1))
      else
         call kdwhereis(z, atmos%nodes, atmos%z, 1, i)
         ans = atmos%T(i) + atmos%b(i)*(z-atmos%z(i))
      endif
      cvh2temp = ans
      end

!---------------------------------------------
      real*8 function cthick2h(t)
      use modAtmosDef
      implicit none
! #include  "Zatmos.h"
      real*8 t       ! input. air thickness in kg/m^2

      real*8 logt, ans
      integer i
      real*8  dod0, fd, a
      
      logt = log(t)
      if(t .ge. atmos%cumd(1) ) then
         ans = atmos%z(1) -
     *      (logt - atmos%logcumd(1) )*atmos%H(1)
      elseif(t .le. atmos%cumd(atmos%nodes)) then
         ans = atmos%z(atmos%nodes) -
     *       Hinf* log(t/atmos%cumd(atmos%nodes))
      else
         call kdwhereis(t, atmos%nodes,  atmos%cumd,  1, i)
!            i is such that     X(i) > x  >= x(i+1);
!            i=1, 2, .. nodes-1
         dod0 =( atmos%cumd(i) - t )/ atmos%d0(i)
         a = atmos%a(i)
         fd = 1. - dod0
         if(a .ne. 0.) then
            ans = (fd**(-a)- 1.0d0)*atmos%H(i)/a + atmos%z(i)
         else
            ans =  -log(fd)* atmos%H(i) + atmos%z(i)
         endif
!/////////////
!       write(0,*) ' t=',t, ' z=', ans, ' i=',i, ' dod0=', dod0, ' a=',a
!       write(0,*) ' atmos.H,z=', atmos.H(i), atmos.z(i)
!/////////////
      endif

      cthick2h = ans
      end

!---------------------------------------------
      real*8 function cthick2den(t)
      use modAtmosDef
      implicit none
! #include  "Zatmos.h"
      real*8 t       ! input. air thickness in kg/m^2

      real*8 logt, ans, temp
      real*8  cthick2h, cvh2den
      

      logt = log(t)
      if(t .gt. atmos%cumd(1) ) then
         temp = cthick2h(t)
         ans =  cvh2den( temp )
      elseif(t .lt. atmos%cumd(atmos%nodes)) then
         temp = cthick2h(t)
         ans =  cvh2den( temp )
      else
         temp = cthick2h(t)
         ans =  cvh2den( temp )
      endif
      cthick2den = ans
      end
!      -------------------------------------
       real*8 function cvh2denp(z)
!      -------------------------------------
!     d rho/dz
       use modAtmosDef
       implicit none
!  #include  "Zatmos.h"
       real*8 z

       real*8  ans
       integer i
       real*8  a  

       if( z .gt. atmos%z(atmos%nodes) ) then
          ans =  0.
       elseif(z .lt. atmos%z(1)) then
          ans =- atmos%rho(1)/atmos%H(1)*
     *          exp( (atmos%z(1)-z)/atmos%H(1) )
       else
         call kdwhereis(z, atmos%nodes, atmos%z, 1, i)
         a = atmos%a(i)
         if(a .ne. 0.d0) then
           ans = atmos%rho(i)*a/atmos%H(i) *
     *      (1.0d0 + a*(z-atmos%z(i))/atmos%H(i) )**(-2.d0-1.0d0/a)
         else
           ans =- atmos%rho(i)/atmos%H(i)* 
     *        exp(- (z-atmos%z(i)/atmos%H(i)))
         endif
       endif
       cvh2denp = ans
       end
!      ----------------------------------
       real*8 function cvh2scaleh(z)
!     ----------------------------------
       use modAtmosDef
       implicit none
!  #include  "Zatmos.h"
       real*8 z
       real*8 ans
       integer i  		   

       if( z .gt. atmos%z(atmos%nodes-1) ) then
          ans = atmos%H(atmos%nodes-1)
       elseif(z .lt. atmos%z(1)) then
          ans = atmos%H(1)
       else
         call kdwhereis(z, atmos%nodes, atmos%z, 1, i)
         ans = atmos%H(i) + atmos%a(i)*(z-atmos%z(i))
       endif
       cvh2scaleh = ans
       end
!      -------------------------------------
       real*8 function cvh2den2p(z)
!     -------------------------------------
!     d(d rho/dz)/dz
       use modAtmosDef
       implicit none
!  #include  "Zatmos.h"
       real*8 z
       real*8 ans
       integer i
       real*8 a

       if( z .gt. atmos%z(atmos%nodes) ) then
          ans = 0.
       elseif(z .lt. atmos%z(1)) then
          ans = atmos%rho(1)/atmos%H(1)/atmos%H(1)*
     *          exp( (atmos%z(1)-z)/atmos%H(1) )
       else
         call kdwhereis(z, atmos%nodes, atmos%z, 1, i)
         a = atmos%a(i)
         if(a .ne. 0.) then
            ans = atmos%rho(i)* (a/atmos%H(i))**2 *
     *      (1. + a*(z-atmos%z(i))/atmos%H(i))**(-3.-1./a)
	 else
            ans = atmos%rho(i)/atmos%H(i)**2  *
     *         exp(-(z-atmos%z(i))/atmos%H(i))
         endif
       endif
       cvh2den2p = ans
       end
!      ---------------------------------------
       real*8 function cvh2thick(z)
!     ---------------------------------------
       use modAtmosDef
       implicit none
!  #include  "Zatmos.h"
       real*8 z

       real*8 ans
		  
       integer i 
       real*8  a

       if( z .gt. atmos%z(atmos%nodes) ) then
          ans = atmos%cumd(atmos%nodes) *
     *     exp((atmos%z(atmos%nodes) - z)/Hinf )
       elseif(z .lt. atmos%z(1)) then
          ans = atmos%cumd(1)*
     *          exp( (atmos%z(1)-z)/atmos%H(1) )
       else
          call kdwhereis(z, atmos%nodes, atmos%z, 1, i)
          a = atmos%a(i)
          if(a .ne. 0.) then
            ans = atmos%cumd(i) - atmos%d0(i)*(1.-
     *          (1+ a*(z-atmos%z(i))/atmos%H(i))**(-1/a))
          else
            ans = atmos%cumd(i) - atmos%d0(i)*(1.-
     *         exp(-(z-atmos%z(i))/atmos%H(i)))
          endif
       endif
       cvh2thick = ans
       end