!!         integer nxsec
!!         real*8  xsecstep    ! making table with 2 mb step
!!         real*8  xsecmin     ! from xsecmin to step xsecstep
!!         parameter (nxsec=66, xsecstep=2.d0)  ! 46-->66; version 8.77
!!         parameter (xsecmin = 1.d0)   !  10 -> 1; version 8.77

       type element 
       sequence
           real*8  A   ! mass number
           real*8  Z   ! atomic number
!!           real(8):: No
!!           real(8):: OrigNo
!!           real(8):: smuNo
!!           real(8):: nsigma
!!           real(8):: w
!!           real(8):: npercm3
!           real*8  sigma(nxsec)  ! 10 mb to 120 mb pp xsec is converted
!                                 !  into this target cross-section
            integer N   ! nucleon number.  N-Z is number of neutrons.
       end type element
