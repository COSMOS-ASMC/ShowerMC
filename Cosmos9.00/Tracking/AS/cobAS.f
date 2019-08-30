!     ******************************************************************
!     *                                                                *
!     * cobAS: creates air shower from a given electron and observes it*
!     *                                                                *
!     ******************************************************************
!
!
      subroutine cobAS(el)
      implicit none
#include "Ztrack.h"
! #include "Zmagfield.h"      
#include "Zobs.h"
#include "Zobsv.h"
#include "Zelemagp.h"
#include "ZmediaLoft.h"
      type(track)::el    ! an electron to produce A%S

      integer l, j
      integer never
      real*8  zobas, zp,  clenbetween2h, length, clen2thick
      real*8  t, elog, eno, age
      real*8  enoc/1./
! /////////
      real*8 sintemp1, sintemp2
! ///////
      logical downgoing
      integer lb, le, ls
      save never
      data never/0/
      real(8):: X0  ! kg/m2

      zp = el%pos%depth      ! depth of electron
      downgoing = el%vec%coszenith .gt. 0.d0
      if(downgoing) then
         lb = 1
         le = NoOfASSites
         ls = 1
      else
         lb = NoOfASSites
         le = 1
         ls = -1
      endif
      do   l=lb, le, ls
!           vertical depth kg/m2
        zobas=ASObsSites(l)%pos%depth

        X0 =Media(MediaNO)%X0g*10.d0 ! X0g->g/cm2

        if( (zp .lt. zobas  .and. downgoing)
     *   .or. ( zp .gt. zobas .and. .not. downgoing) ) then
!            electron will cross A.S site.
!               thickness to r.l
            if(abs(el%vec%coszenith) .gt. 0.5) then
               t=abs(zobas-zp)/X0/abs(el%vec%coszenith)
            else
!                  error check: since A.S plane is
!                  spherical surface, particle may not pass
!                  the next surface.  (if completely horizontal,
!                  it should pass)
               sintemp1 = sqrt(1.d0-  el%vec%coszenith**2)
               sintemp2 = el%pos%radiallen * sintemp1/
     *                    ASObsSites(l)%pos%radiallen
               if( sintemp2 .le. 1.d0) then
!                 slant length between two depths
                  length = clenbetween2h(el%pos%radiallen, 
     *                  ASObsSites(l)%pos%radiallen,
     *                  el%vec%coszenith)
                  t = clen2thick(el%pos%height, el%vec%coszenith,
     *                        length) / X0
               else
!                   particle dose not cross the surface. neglect A.S
                  CompASNe(l) = 0.
                  CompASAge(l) = 0.
                  goto 500
               endif
            endif
!                log of electron energy in unit of critical energy
            elog = log10(el%p%fm%p(4)/Media(MediaNo)%Ecrit) 
!                        !-get ne and age
            call cNeByApproxB(3, elog, t,   eno, age)
!                 count ne
!                 if thins=T, wgt may be diff. from 1.
            eno=eno*el%wgt
!
            CompASNe(l) = eno
            CompASAge(l) = age
!
            ASObsSites(l)%esize=ASObsSites(l)%esize + eno
            ASObsSites(l)%age = ASObsSites(l)%age + eno*age
            if(eno .lt. enoc) then
!             ne is very small
               if( age .gt. 1.) then
!                     old so that cannot become larger; discard it
                  do j = l+ls, le, ls
                     CompASNe(j) = 0.
                     CompASAge(j) = 2.
                  enddo
                  goto 600
               else
!                    age is young, a.s might become larger
               endif
            else
!                  ne is large.  if needed, lateral obs.
            endif
         else
            CompASNe(l) = 0.
            CompASAge(l) = 0.
         endif
 500     continue
      enddo
  600 continue
      if(never .ne. 1) then
         call chookHybAS(el, never)
      endif
      end
