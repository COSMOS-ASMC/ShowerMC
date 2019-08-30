      module BPKochMotz
      implicit none
!      private
      real(8)::  Z13, Z23, ame, apme, delta, cf,
     *        lnZ43,  lnZ83, Z, Z2, BHnorm, NonScEme,
     *        al183z, ccz, bcoef, x0g, fz
!       x0g:  r.l for virtual matte:r A=1 and Z.
!          next is exception.  dependent on Eg or Ee
      real(8):: Ebyme    ! Ee/Me or  Eg/Me


      contains
      function epBrem() result(ans)
 !     use  BPKochMotz
      implicit none
      real(8)::ans
!         
      real(8)::  bigf
      ans =  bigf()*cf
      end function epBrem

!
      function bigf( ) result(ans)
!      use  BPKochMotz
      implicit none
      real(8):: ans
      ans = 1000.* Z2
      end function bigf

      end module BPKochMotz
