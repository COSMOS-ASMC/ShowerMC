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
      
      write(io, '(1p, 6g13.4)'
     *  aSite%pos%depth*0.1), aSite%pos%height, 
     *  aSite%pos%radiallen/1000., 
     * aSite%pos%xyz%r(1:3)
      end
!    *********************************print header
      subroutine cprObsSiteHd(io)
!    **************************************
      integer io ! inpu%t dev. no.
      write(io, '(a)')
     * ' depth (g/cm2)  Height(m)  R(km) and'//
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
      
      write(io, '(1p,7g13.4)')
     *  aSite%pos%depth*0.1, aSite%pos%height, 
     *  aSite%pos%radiallen/1000.,aSite%mu, 
     *  aSite%pos%xyz%r(1:3)
      end
!    *********************************print header
      subroutine cprASObsSiteHd(io)
!    **************************************
      integer io ! input.
      write(io, '(a)' ) 
     * ' depth (g/cm2)    Height(m)   Distance to E-center(km)'//
     * ' Molere Unit(m)     x,y,z in XYZ system (m)'
      end



