!        pi-n total xsection , elastic xs
      subroutine cpimnTotXs(p, xs)
      implicit none
      real*8 p  ! input.  momentum of n. in GeV
      real*8 xs     ! output. total np cross section in mb
      call cpippTotXs(p, xs)
      end      subroutine cpimnTotXs
!         
      subroutine cpimnElaXs(p, xs)
!           pi+ p elastic cross section in mb
      implicit none
      real*8 p ! input.  momentum of n in GeV
      real*8  xs   ! output np elastic xs. mb.
      call cpippElaXs(p, xs)
      end      subroutine cpimnElaXs
      subroutine cpimnInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs

      real(8)::txs, exs
      call cpippInelaXs(p, xs)
      end      subroutine cpimnInelaXs

