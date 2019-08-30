!  This is to express atmosphere using cspline function 
!
!      -------------------------------------
       real*8 function cvh2den(z)
!     --------------------------vertical height to density
       use modAtmosDef
       implicit none
!  #include  "Zatmos.h"
       real*8 z  ! input. m
!       real*8 zsave
       real*8 temp, ans
!       data zsave/-1.d30/
!       save zsave, ans

!       if(z .ne. zsave) then
          if( z .gt. atmos%z(atmos%nodes) ) then
             ans = atmos%rho(atmos%nodes)*
     *            exp(-(z-atmos%z(atmos%nodes))/Hinf)
          elseif(z .lt. atmos%z(1)) then
             ans = atmos%rho(1)*
     *          exp( (atmos%z(1)-z)/atmos%H(1) )
          else
             call kcsplIntp(atmos%z, atmos%logrho, atmos%nodes,
     *         atmos%coefh2r, maxnodes, z, temp)               
             ans = exp(temp)
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

      if( z .gt. atmos%z(atmos%nodes) ) then
         ans = atmos%T(atmos%nodes)
      elseif(z .lt. atmos%z(1)) then
         ans = atmos%T(1) + atmos%b(1)*(z - atmos%z(1))
      else
         call kcsplIntp(atmos%z, atmos%T, atmos%nodes,
     *        atmos%coefh2T, maxnodes, z, ans)
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
      

      logt = log(t)
      if(t .gt. atmos%cumd(1) ) then
         ans = atmos%z(1) -
     *      (logt - atmos%logcumd(1) )*atmos%H(1)
      elseif(t .lt. atmos%cumd(atmos%nodes)) then
         ans = atmos%z(atmos%nodes) -
     *       Hinf* log(t/atmos%cumd(atmos%nodes))
      else
         call kcsplIntp(atmos%logcumdi, atmos%zi, atmos%nodes,
     *      atmos%coefd2h, maxnodes, logt, ans)
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
         call kcsplIntp(atmos%logcumdi, atmos%logrhoi, 
     *      atmos%nodes,
     *      atmos%coefd2r, maxnodes, logt, ans)
         ans = exp(ans)
      endif
      cthick2den = ans
      end
!      -------------------------------------
      real*8 function cvh2denp(z)
      use modAtmosDef
!      -------------------------------------
!          d rho/dz
       implicit none
!  #include  "Zatmos.h"
       real*8 z

       real*8 d2, ans, temp

       if( z .gt. atmos%z(atmos%nodes) ) then
          ans =  0.
       elseif(z .lt. atmos%z(1)) then
          ans =- atmos%rho(1)/atmos%H(1)*
     *          exp( (atmos%z(1)-z)/atmos%H(1) )
       else
          call kcsplDif(atmos%z, atmos%nodes,
     *       atmos%coefh2r, maxnodes,  z,  ans, d2)
          call kcsplIntp(atmos%z, atmos%logrho, atmos%nodes,
     *    atmos%coefh2r, maxnodes, z, temp)               
          ans = exp(temp) * ans
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

       if( z .gt. atmos%z(atmos%nodes-1) ) then
          ans = atmos%H(atmos%nodes-1)
       elseif(z .lt. atmos%z(1)) then
          ans = atmos%H(1)
       else
          call kcsplIntp(atmos%z, atmos%H, atmos%nodes-1,
     *       atmos%coefh2H, maxnodes, z, ans)    
       endif
       cvh2scaleh = ans
       end
!      -------------------------------------
       real*8 function cvh2den2p(z)
!      -------------------------------------
!     d(d rho/dz)/dz
       use modAtmosDef
       implicit none
! #include  "Zatmos.h"
       real*8 z
       real*8 ans, d1, d2, temp

       if( z .gt. atmos%z(atmos%nodes) ) then
          ans = 0.
       elseif(z .lt. atmos%z(1)) then
          ans = atmos%rho(1)/atmos%H(1)/atmos%H(1)*
     *          exp( (atmos%z(1)-z)/atmos%H(1) )
       else
         call kcsplDif(atmos%z, atmos%nodes,
     *       atmos%coefh2r, maxnodes,  z,  d1, d2)
         call kcsplIntp(atmos%z, atmos%logrho, atmos%nodes,
     *    atmos%coefh2r, maxnodes, z, temp)               
          ans = exp(temp) * (d2 + d1**2)
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

       real*8 ans, temp

       if( z .gt. atmos%z(atmos%nodes) ) then
          ans = atmos%cumd(atmos%nodes) *
     *     exp((atmos%z(atmos%nodes) - z)/Hinf )
       elseif(z .lt. atmos%z(1)) then
          ans = atmos%cumd(1)*
     *          exp( (atmos%z(1)-z)/atmos%H(1) )
       else
          call kcsplIntp(atmos%z, atmos%logcumd, atmos%nodes,
     *      atmos%coefh2d, maxnodes, z, temp)
          ans = exp(temp)
       endif
       cvh2thick = ans
       end
      
