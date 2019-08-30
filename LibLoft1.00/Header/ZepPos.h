!             position vector
       type epPos 
       sequence
!       union
!       map
        real*8 x, y, z  !  x,y,z in cm
!       endmap
! if union/map is used to be able to use r(:)
!       assignment like p(n)=epPos(1.0d, 2.d0, 3.d0)
!      becmoes error
!       map
!        real(8)::r(3)
!       endmap
!      end union
       end type epPos
