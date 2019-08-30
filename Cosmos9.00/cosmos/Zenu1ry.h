!        structure need to compute geomagntic effect and
!        transformation between 'enu' and '1ry' systems.
!
         type  enu1ry
           sequence
              real *8 mat(3, 3)     ! transformation matrix 'enu'<->'1ry'
              real *8 b             ! absolute mag. of mag. field.
              real *8 bt            ! b * sin(alfa). 
!                                     alfa = angle between  1ry direction
!                                     and the magnetic field direction.
              real *8 cosa          ! cos(alfa)
              real *8 sina          ! sin(alfa)
         end type enu1ry

