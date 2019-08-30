      subroutine cdet2xyzFE(aTrack, outr, outdir)
!      front-end routine for  converting position and direction cos
!      of aTrack in the detector system into the  E-xyz system.
! Similar routines follow:
!           det->xyz xyz->det
!           prim->xyz, xyz->prim
!           det->prim, prim->det 
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
      type(track)::aTrack  ! input Track info
      real(8),intent(out):: outr(3), outdir(3)

      call cdet2xyz(ObsSites(aTrack%where)%pos%xyz,
     *  aTrack%pos%xyz%r, outr)
      call cdet2xyzD(aTrack%vec%w%r, outdir)
      end

      subroutine cxyz2detFE(aTrack, outr, outdir)
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
      type(track)::aTrack  ! input Track info
      real(8),intent(out):: outr(3), outdir(3)

      call cxyz2det(ObsSites(aTrack%where)%pos%xyz,
     *  aTrack%pos%xyz%r, outr)
      call cxyz2detD(aTrack%vec%w%r, outdir)
      end

      subroutine cprim2xyzFE(aTrack, outr, outdir)
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
      type(track)::aTrack  ! input Track info
      real(8),intent(out):: outr(3), outdir(3)

      call cprim2xyz(ObsSites(aTrack%where)%pos%xyz,
     *  aTrack%pos%xyz%r, outr)
      call cprim2xyzD(aTrack%vec%w%r, outdir)
      end

      subroutine cxyz2primFE(aTrack, outr, outdir)
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
      type(track)::aTrack  ! input Track info
      real(8),intent(out):: outr(3), outdir(3)

      call cprim2xyz(ObsSites(aTrack%where)%pos%xyz,
     *  aTrack%pos%xyz%r, outr)
      call cprim2xyzD(aTrack%vec%w, outdir)
      end


      subroutine cprim2detFE(aTrack, outr, outdir)
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
      type(track)::aTrack  ! input Track info
      real(8),intent(out):: outr(3), outdir(3)

      call cprim2xyz(ObsSites(aTrack%where)%pos%xyz,
     *  aTrack%pos%xyz%r, outr)
      call cxyz2det(ObsSites(aTrack%where)%pos%xyz,
     *      outr, outr)
      call cprim2xyzD(aTrack%vec%w%r, outdir)
      call cxyz2detD(outdir, outdir)
      end

      subroutine cdet2primFE(aTrack, outr, outdir)
      implicit none
#include  "Zglobalc.h"
#include  "Ztrack.h"
! #include  "Zmagfield.h"      
#include  "Ztrackp.h"
#include  "Ztrackv.h"
#include  "Zobs.h"
#include  "Zobsp.h"
#include  "Zobsv.h"
      type(track)::aTrack  ! input Track info in det sys.
      real(8),intent(out):: outr(3), outdir(3)
      type(track)::tTrack  !  temp track

      tTrack=aTrack  ! det
      call cdet2xyz(ObsSites(tTrack%where)%pos%xyz,
     *  tTrack%pos%xyz%r, outr)   ! xyz
      call cxyz2prim(ObsSites(tTrack%where)%pos%xyz,
     *   outr, outr)  ! prim
      call cdet2xyzD( aTrack%vec%w%r, outdir)
      call cxyz2primD(outdir, outdir)
      end

