      implicit none
!      
!       read skelton data and show it;
! This is  another version of showing skelton data in ascii
! (see seeascii.f)
!
#include  "Ztrack.h"
      include "Zprivate.h"
      type(ob)::oo
      type(child)::cc
      type(parent)::pp 
      integer i, nlow, cumnum, num, ir(2)
      type(track)::zf

      open(21, file='localhost_skelnode', form='unformatted')
      do while(.true.)
         read(21, end=100) cumnum, num, ir, zf
         write(*,'(a, 4i11)' ) 'h ', cumnum, num, ir
         write(*,'(a, i3,2i5,4x,g12.3,3x,f8.5,g12.3)')
     *    'f ',
     *    zf%p%code, zf%p%subcode, zf%p%charge,
     *    zf%p%fm%p(4), zf%vec%coszenith,
     *    zf%pos%depth/10.

         read(21) Np
         write(*,'(a, i6)' ) 'nob ', Np
         do i = 1, Np
            read(21) oo

            write(*,'(4i3,9g14.6)') oo%where, oo%code, oo%subcode, 
     *        oo%charge, oo%atime, oo%erg, 
     *        oo%mass, oo%x, oo%y, oo%wx, oo%wy, oo%wz, 
     *        oo%zenith
         enddo
         nlow =1
         do while (nlow .gt. 0)
            read(21) nlow, pp
            if(nlow .gt. 0) then
               write(*,'(a,i8)') 'nl ', nlow
               if(nlow .gt. 0) then
                  write(*,'(a, 8g13.6,3i3,g13.6)')
     *                 'p ', pp%posx, pp%posy, pp%posz,pp%coszenith,
     *                 pp%depth, pp%colHeight, 
     *                 pp%height, pp%atime, pp%where, pp%code,
     *                 pp%asflag, pp%erg
                  do i = 1, nlow
                     read(21) cc
                     write(*,'(a, 3i3,5g13.6)')
     *                    'c ', cc%code, cc%subcode, cc%charge, 
     *                    cc%fm(1), cc%fm(2),cc%fm(3), cc%fm(4),
     *                    cc%mass 
                  enddo
               endif
            endif
         enddo
      enddo
 100  continue
      end
