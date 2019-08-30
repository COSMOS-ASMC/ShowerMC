  module modDCS
    use modTPXS
    implicit none
!    integer,parameter::nmu=111  ! thinned version
    integer,parameter::nmu=606   ! non thinned version
    real(8)::muval(nmu) = -1.0
    real(8)::logmuval(nmu)   ! logmuval(1) is log( muval(2)/5)
    type DCSconst
!     This is used only when making sampling table or to read
!         existing data.
       real(8):: dcs(nmu,nEneg)  !  dcs(idxmu, idxE)
                              ! size is dcs(606, 106).  For 20 media, ~ 606x106*8*20= 10 MB
       real(8):: logdcs(nmu, nEneg)
!       cubic spline coef for mu (@fixed grid Energy)
!!cs       real(8):: coefdcs(nmu-1, 3, nEneg)
!  For these 3, (20+60)*2 = 160 MB
    end type DCSconst
  end module modDCS
