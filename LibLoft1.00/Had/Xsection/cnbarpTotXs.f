!        nbarp total xsection, nbarp ela
      subroutine cnbarpTotXs(p, xs)
      implicit none
      real*8 p  ! input.  momentum of nbar. in GeV
      real*8 xs     ! output. total nbarp cross section in mb
      call cpbarpTotXs(p, xs)
      end      subroutine cnbarpTotXs
      subroutine cnbarpElaXs(p, xs)
!           pbarp elastic cross section in mb
      implicit none
      real*8 p ! input.  momentum of pbar in GeV
      real*8  xs   ! output pbarp elastic xs. mb.
      call cpbarpElaXs(p, xs)
      end      subroutine cnbarpElaXs
      subroutine cnbarpInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs

      real(8)::txs, exs
      call cnbarpTotXs(p, txs)
      call cnbarpElaXs(p, exs)
      xs =max( txs - exs, 0.d0)
      end      subroutine cnbarpInelaXs

