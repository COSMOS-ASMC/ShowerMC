      subroutine epLightPushCells(io, fmt)
      use modepLightPty
      use modepLightEdepo
      use modepLight
      implicit none
#include  "Zep3Vec.h"
#include  "ZepTrack.h"
#include  "Zcnfig.h"
#include  "Zcode.h"
!  #include  "ZepManager.h"

      integer,intent(in)::io  !  if 0, push cells into stack else write
                                !  it to file# "io"
      integer,intent(in)::fmt   ! user when io /= 0.     >0  --> formatted
                                !    < 0  --> unformatted
      integer::i, j
      integer::clcmp
      integer::  l, m, n, lm, mm, nm

      real(8), pointer::dx, dy, dz
      integer, pointer::nx, ny, nz

       type(epTrack)::  cell

      real(8)::epCellIdx2v

        ! function to convert cell index to boundary value

      do i = 1, LightCompNo
         clcmp = i
         if(  allocated( parts(i)%edepo ) ) then
            j =Lcomp(i)%compno
!/////////////
            if(j <=  0 ) then
               write(0,*) ' cn=', j, ' for lightcomp #=',i
            endif
!///////////////////

            dx => Lcomp(i)%dx
            dy => Lcomp(i)%dy 
            dz => Lcomp(i)%dz
            nx => Lcomp(i)%nx
            ny => Lcomp(i)%ny
            nz => Lcomp(i)%nz

            if( Det%cmp(j)%struc(1:3) == "box" .or. 
     *          Det%cmp(j)%struc(1:7) == "octagon"  ) then
               lm = 0
               mm = 0
               nm = 0
            elseif( Det%cmp(j)%struc(1:3) == "cyl" .or.
     *              Det%cmp(j)%struc(1:4) == "ecyl" ) then
               lm = -nx
               mm = -ny
               nm = 0
            elseif( Det%cmp(j)%struc(1:4) == "pipe"  ) then
               lm = -nx
               mm = -ny
               nm = 0
            else
               write(0,*)  Det%cmp(j)%struc,
     *                 ' is not yet supported for Light'
               stop
            endif
            do l = lm, nx
               do m = mm, ny
                  do n = nm, nz
                     cell%p%fm%p(4) =parts(clcmp)%edepo(l, m, n) 
                     if( cell%p%fm%p(4) > 0.) then
                        cell%cn = j
                        cell%pos%x = epCellIdx2v(l,dx)
                        cell%pos%y = epCellIdx2v(m,dy)
                        cell%pos%z = epCellIdx2v(n,dz)
                        cell%p%code = kEdepo
                        cell%p%subcode = 0
                        cell%p%charge = 0
                        cell%w%x = dx
                        cell%w%y = dy
                        cell%w%z = dz
                        cell%wl = 0. ! dummy in this case
                        cell%p%mass = 0. ! dummy in this case
                        if(io == 0) then
                           call eppush(cell)
                        else
       !   common to Edepo and chgPath
       !           code, charge, E, pos, dir,  mass, Cn 
       !   for Edepo:    charge=0. dir=(dx,dy,dz)= cell size
       !   for chgPath:            dir=dl(dx,dy,dz)= directed path length
       !   +primary file:
       !    mulTExyzdirmassCn  (xyz means  pos is local coord since Cn is given)
                           if(fmt > 0 ) then
                              write(io, 
     *          '(i6,2i4,1p,4g15.7,3g18.9,2g12.4,0p,i6)' )
     *                        kEdepo, 
     *                        cell%p%subcode,
     *                        cell%p%charge,
     *                        cell%p%fm%p(4), cell%pos, cell%w,
     *                        cell%wl, cell%p%mass,  j
                           else
                              write(io)  kEdepo, 
     *                        int(cell%p%subcode),
     *                        int(cell%p%charge),
     *                        cell%p%fm%p(4), cell%pos, cell%w,
     *                        cell%wl, cell%p%mass,  j
                           endif
                        endif
                     endif
                  enddo
               enddo
            enddo
            deallocate( parts(i)%edepo ) 
         endif
      enddo

      end  subroutine epLightPushCells
