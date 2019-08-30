      subroutine cfxTgtChg(ia, iz, tcg)
!           fix target charge (for nucleon target inside 'A')
      implicit none
      integer ia, iz    !   input target mass and charge no.
      integer tcg       !   ouput. selected nucleon's charge

      real*8 u

      call rndc(u)
      if(u .le.  float(iz)/ia) then
         tcg=1
      else
         tcg=0
      endif
      end
