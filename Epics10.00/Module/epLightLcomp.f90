!     not used ???????????
module modepLightLcomp
  use modepLightPty
  !        for each  light component, we define.
  type(compPty)::Lcomp(maxUniqLightCMP)
!ontains
  subroutine epLightManip(i)  
    implicit none
    integer,intent(in)::i
    Lcomp(i)%refracIndex = 10.
  end subroutine epLightManip
end module modepLightLcomp


