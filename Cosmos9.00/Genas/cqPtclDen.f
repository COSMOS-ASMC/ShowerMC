      subroutine cqPtclDen(lat, depindx, how, r, rho)

!      Since the observed particle density is heavily affected by
!      the detector structure, we assume here some detecotrs for 
!      certain input parameters.

      implicit none
      integer lat   ! input. 1 --> Lateral distribution used in GENAS
                    !              is employed. (recomended)
                    !        2 --> nkg function is used, making the
                    !              Moliere unit length half.
                    !        3 --> nkg function is used with normal
                    !              Moriele unit
                    !        4 --> Bare electron lateral distribution
                    !              is used which created by an electron
                    !              primary.   
      integer depindx ! input. A%S observation depth index.
                      !        1 to NoOfASSites.

      integer how   ! input.  has sence, if lat = 1.
                    !        1--> very normal scintillator detector is
                    !             assumed; thin iron + scintillator
                    !             of  ~ 4cm thick. 
                    !             
                    !        2--> 0.5 cm lead plate is added to the
                    !             how=1 case.
                    !   
      real*8 r      ! input. distance from the core in m (pependicular
                    !        to the axis) where you want to have the
                    !        effective particle density.

      real*8 rho    ! output. effective average particle number in an area of
                    !         1 m2 at r.  The area is assumed to be
                    !         horiazontal (if ObsPlane = 1) or perpendicular
                    !         to the shower direction (if ObsPlane =2)

#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
#include "Zcode.h"

      type(track):: inci
      
      type(coord):: angle
      real*4  s, e0, cosz, rmu, rhog
      real*8  hr, cnkg, mu, s8

      call cqIncident(inci, angle)
      s = ASObsSites(depindx)%age
      e0 = inci%p%fm%p(4) /1000.   ! TeV
      cosz = inci%vec%coszenith  
      mu = ASObsSites(depindx)%mu   
      rmu = r /mu   ! in Moriele unit 

      if(lat .eq. 1) then
         if(inci%p%code .eq. kphoton .or. 
     *      inci%p%code .eq. kelec) then
            if(how .eq. 1) then
               call kdig0(e0, cosz, s, rmu, rhog)
            elseif(how .eq. 2) then
               call kdigb0(e0, cosz, s, rmu, rhog)
            else
               call cerrorMsg('how is wrong in cqPtclDen.f', 0)
            endif
         else
            if(how .eq. 1) then
               call kdip0(e0, cosz, s, rmu, rhog)
            elseif(how .eq. 2) then
               call kdipb0(e0, cosz, s, rmu, rhog)
            else
               call cerrorMsg('how is wrong in cqPtclDen.f', 0)
            endif
         endif
      elseif(lat .eq. 2) then
         hr = rmu*2.
         s8 = s
         rhog = cnkg(s8, hr)
      elseif(lat .eq. 3) then
         hr = rmu
         s8 = s
         rhog = cnkg(s8, hr)
      elseif(lat .eq. 4) then
         call klee(s, rmu, rhog)
      else
         call cerrorMsg('lat is wwong in cqPtclDen', 0)
      endif
      rho = rhog * mu * mu
      if(ObsPlane .eq. 1) then
         rho =rho * cosz
      endif
      end

         




