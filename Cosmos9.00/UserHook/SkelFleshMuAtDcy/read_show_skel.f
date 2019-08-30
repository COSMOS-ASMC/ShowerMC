      implicit none
!      
!       read skelton data and show it;
! This is  another version of showing skelton data in ascii
! (see seeascii.f)
!
#include  "Ztrack.h"
      include "Zprivate.h"
      type(ob):: oo
      type(child):: cc
      type(parent):: pp 
      integer i, nlow, cumnum, num, ir(2)
      type(track)::zf
      real*8 user
      real*4 userE(2)
      integer*2 userI(4)
      equivalence (user, userE(1)), (user, userI(1))
      logical kbitest
      integer  codex

      open(21, file='localhost_skelnode', form='unformatted')
      do while(.true.)
         read(21, end=100) cumnum, num, ir, zf
         write(*,'(a, 4i11)' ) 'h ', cumnum, num, ir
         write(*,'(a, i3,2i5,4x,g12.3,3x,f8.5,g12.3)')
     *    'f ',
     *    zf.p.code, zf.p.subcode, zf.p.charge,
     *    zf.p.fm.p(4), zf.vec.coszenith,
     *    zf.pos.depth/10.

         read(21) Np
         write(*,'(a, i6)' ) 'nob ', Np
         do i = 1, Np
            read(21) oo
            user=oo.user
            if( kbitest(userE(2), 1) ) then
               codex = 4
            else
               codex = 6
            endif
            write(*,'(g14.6, i3,3g11.3)') 
     *         oo.erg, codex, userI(1)/1000.,  userI(2)/1000.,
     *         userE(2)
!            write(*,
!        *      '(4i3, 1p4E11.3, 0p, 2f8.4,f10.6, i3, 1p g11.3)')
!        *        ldep,  code,  aTrack.p.subcode, aTrack.p.charge,
!        *        Ek, aTrack.t,
!        *        aTrack.pos.xyz.x, aTrack.pos.xyz.y,
!        *        -aTrack.vec.w.r(1),  -aTrack.vec.w.r(2),  wz,
!        *        userI(1), userE(2)
                                          
         enddo
         nlow =1
         do while (nlow .gt. 0)
            read(21) nlow, pp
            if(nlow .gt. 0) then
               write(*,'(a,i8)') 'nl ', nlow
               if(nlow .gt. 0) then
                  user = pp.user
                  write(*,'(a, 8g13.5,3i3,4g13.6)')
     *                 'p ', pp.posx, pp.posy, pp.posz,pp.coszenith,
     *                 pp.depth, pp.colHeight, 
     *                 pp.height, pp.atime, pp.where, pp.code,
     *                 pp.asflag, pp.erg, userI(1)/1000., 
     *                 userI(2)/1000., userE(2)
                  do i = 1, nlow
                     read(21) cc
                     user = cc.user
                     write(*,'(a, 3i3,5g13.6, 3g11.3)')
     *                    'c ', cc.code, cc.subcode, cc.charge, 
     *                    cc.fm(1), cc.fm(2),cc.fm(3), cc.fm(4),
     *                    cc.mass, userI(1)/1000., userI(2)/1000.,
     *                    userE(2)
                  enddo
               endif
            endif
         enddo
      enddo
 100  continue
      end

