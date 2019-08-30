!        pi+n total xsection , elastic xs
      subroutine cpipnTotXs(p, xs)
      implicit none
      real*8 p  ! input.  momentum of n. in GeV
      real*8 xs     ! output. total np cross section in mb
      call cpimpTotXs(p, xs)
      end      subroutine cpipnTotXs
!         
      subroutine cpipnElaXs(p, xs)
!           pi+ p elastic cross section in mb
      implicit none
      real*8 p ! input.  momentum of n in GeV
      real*8  xs   ! output np elastic xs. mb.
      call cpimpElaXs(p, xs)
      end      subroutine cpipnElaXs
      subroutine cpipnInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs
      call cpimpInelaXs(p, xs)
      end subroutine cpipnInelaXs

