!          print an observation site information
!
      subroutine cprObsSite(io, aSite)
      implicit none
#include "Zmagfield.h"
#include "Zcoord.h"
#include "Zdirec.h"
#include "Zpos.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
      integer io ! input. dev. no.
      type(site)::aSite
      
      write(io, *)sngl(aSite%pos%depth*0.1),
     * sngl(aSite%pos%height), 
     * sngl(aSite%pos%radiallen/1000.),
     * aSite%pos%xyz%r(1), aSite%pos%xyz%r(2), aSite%pos%xyz%r(3)
      end
!    *********************************print header
      subroutine cprObsSiteHd(io)
!    **************************************
      integer io ! inpu%t dev. no.
      write(io, *)
     * ' depth (gr/cm2)  Height(m)  Distance to E-center(km)'//
     * '   x,y,z in XYZ system(m)'
      end
! 
!          print an observation site information
!
      subroutine cprASObsSite(io, aSite)
      implicit none
#include "Zmagfield.h"
#include "Zcoord.h"
#include "Zdirec.h"
#include "Zpos.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
      integer io ! input . output dev. no.
      type(assite)::aSite
      
      write(io, *)sngl(aSite%pos%depth*0.1),
     * sngl(aSite%pos%height), 
     * sngl(aSite%pos%radiallen/1000.),
     * sngl(aSite%mu), 
     * aSite%pos%xyz%r(1), aSite%pos%xyz%r(2), aSite%pos%xyz%r(3)
      end
!    *********************************print header
      subroutine cprASObsSiteHd(io)
!    **************************************
      integer io ! input.
      write(io, *)
     * ' depth (gr/cm2)    Height(m)   Distance to E-center(km)'//
     * ' Molere Unit(m)     x,y,z in XYZ system (m)'
      end



