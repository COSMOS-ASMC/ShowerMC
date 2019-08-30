      subroutine ckmnTotXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs
      call ckmpTotXs(p,xs)
      end      subroutine ckmnTotXs
      subroutine ckmnElaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs
      call ckmpElaXs(p,xs)
      end      subroutine ckmnElaXs
      subroutine ckmnInelaXs(p, xs)
      implicit none
      real(8),intent(in)::p
      real(8),intent(out)::xs

      real(8)::txs, exs
      call ckmnTotXs(p,txs)
      call ckmnElaXs(p, exs)
      xs =max( txs - exs, 0.d0)
      end      subroutine ckmnInelaXs

