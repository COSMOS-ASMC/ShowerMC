!      real*8 x, y, z, t, f
!            real*8 pi,   
!     *             Todeg 
!        parameter(pi = 3.141592653589793238d0, 
!     *    Todeg = 180.d0/pi)c

!        do while(.true.)
!           read(*, *) t, f
!          if(f .gt. 500.) stop
!           x = sin(t/Todeg) *  cos(f/Todeg)
!           y = sin(t/Todeg) *  sin(f/Todeg)
!           z = cos(t/Todeg)
!           call kdir2deg(x, y, z, t, f)
!           write(*, *) ' t, f=', t, f
!        enddo
!      end
        
      subroutine kdir2deg(dx, dy, dz, theta, fai)
      implicit none
      real*8 dx  ! input direction cos x comp.
      real*8 dy  ! input //            y comp.
      real*8 dz  ! input //            z comp.

      real*8 theta  ! output.  zenith angle in deg. 0 to 180
      real*8 fai    ! output.  azimuthal angle in deg. 0 to 360

!              constants thru Cosmos
            real*8 pi,   
     *             Todeg 
        parameter(pi = 3.141592653589793238d0, 
     *    Todeg = 180.d0/pi)

        theta = acos(dz)*Todeg
        if(dx .eq. 0) then
           fai = 0.
        else
           fai = atan2(dy, dx)*Todeg
           if(fai .lt. 0.) then
              fai = fai + 360.
           endif
        endif
      end
