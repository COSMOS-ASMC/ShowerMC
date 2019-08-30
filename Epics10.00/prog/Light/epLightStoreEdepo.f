!      when Light=21 or 12, energy deposit is stored in a
!      number of cells in each "light" component.
!      This is used to store the deposit in the cell.
      function epCellIdx(v, dv)
      integer::epCellIdx

      real(8),intent(in)::v
!      real(8),pointer,intent(in)::dv
      real(8),intent(in)::dv
      if( v >=0. ) then
         epCellIdx = ceiling(v/dv)
      else
         epCellIdx = floor(v/dv)
      endif
      end

      function epCellIdx2v(idx, dv)
 !           cell idx  to boundary value  conversion 
      real(8)::epCellIdx2v
      integer,intent(in)::idx
!      real(8),pointer,intent(in)::dv  ! digitizer unit
      real(8),intent(in)::dv  ! digitizer unit
      if(idx > 0 ) then
          ! idx  1       2       3   
          !  (0, dx], (dx,2dx], (2dx,3dx]
         epCellIdx2v = (idx - 1)* dv
      elseif(idx < 0 ) then
          ! idx  -1       -2       3   
          !  [-dx, 0), [-2dx,-dx) [-3dx,-2dx)
         epCellIdx2v = (idx + 1) * dv
      else
         epCellIdx2v = 0.
           !  idx=0  <--> 0.
      endif
      end

      subroutine epLightStoreEdepo(info, aTrack)
      use modepLight
      use modepLightPty
      use modepLightEdepo
        !
      implicit none
#include "ZepTrackp.h"
#include "ZepTrackv.h"
#include "Zcnfig.h"
#include "ZsepManager.h"
       !  
      integer,intent(in):: info ! 0--> track with Move.dl length
                                ! 1--> dl is almost 0
       type(epTrack)::  aTrack  ! input. track bef. move

      real(8), pointer::dx, dy, dz
      integer, pointer::nx, ny, nz
      real(8):: x, y, z

      integer:: epCellIdx    ! function to get cell index

      integer::nnx, nny, nnz 
      real(8)::l, cl, hl, dl, dE      
!////////
!      call Lcompchk(" Zb", cLcompNo) 
!////////
      dy => Lcomp(cLcompNo)%dy
      dx => Lcomp(cLcompNo)%dx
      dz => Lcomp(cLcompNo)%dz
      nx => Lcomp(cLcompNo)%nx
      ny => Lcomp(cLcompNo)%ny
      nz => Lcomp(cLcompNo)%nz

      x = aTrack%pos%x
      y = aTrack%pos%y
      z = aTrack%pos%z

      if( .not.  allocated( parts(cLcompNo)%edepo )) then
         if( Det%cmp(Cn)%struc(1:3) == "box" ) then
            allocate( parts(cLcompNo)%edepo(0:nx,0:ny,0:nz) ) 
         elseif( Det%cmp(Cn)%struc(1:7) == "octagon" ) then
            allocate( parts(cLcompNo)%edepo(0:nx,0:ny,0:nz) ) 
 !             index 0 will not be used for cyl /ecyl 
         elseif( Det%cmp(Cn)%struc(1:3) == "cyl" ) then
            allocate( parts(cLcompNo)%edepo(-nx:nx, -ny:ny, 0:nz))

         elseif( Det%cmp(Cn)%struc(1:4) == "ecyl" ) then
            allocate( parts(cLcompNo)%edepo(-nx:nx, -ny:ny, 0:nz))

         elseif( Det%cmp(Cn)%struc(1:4) == "pipe" ) then
            allocate( parts(cLcompNo)%edepo(-nx:nx, -ny:ny, 0:nz))

         else
            write(0,*) ' struc=',Det%cmp(Cn)%struc
            write(0,*) 'not supported for epLightEdepo'
            stop
         endif
         parts(cLcompNo)%edepo(:,:,:) = 0.
      endif
      
      if( info  ==  1 ) then
!          very short point 
         nnx = epCellIdx(x, dx)
         nny = epCellIdx(y, dy)
         nnz = epCellIdx(z, dz)
!         if(nnx .ge. 1 .and. nnx .le. nx ) then
!            if(nny .ge. 1 .and. nny .le. ny ) then
!               if(nnz .ge. 1 .and. nnz .le. nz ) then
!                  keep in GeV
         parts(cLcompNo)%edepo(nnx, nny, nnz) =
     *        parts(cLcompNo)%edepo(nnx, nny, nnz) + Move%dEeff
!               endif
!            endif
!         endif
      else
         l = 0.
         dl = Lcomp(cLcompNo)%dmin 
         hl = dl/2
         dE = Move%dEeff * dl / Move%dl
         do while ( l + dl  < Move%dl )
            cl = l + hl
            x = aTrack%pos%x + cl*aTrack%w%x 
            nnx = epCellIdx(x, dx)
            y = aTrack%pos%y + cl*aTrack%w%y 
            nny = epCellIdx(y, dy)
            z = aTrack%pos%z + cl*aTrack%w%z
            nnz = epCellIdx(z, dz)

            parts(cLcompNo)%edepo(nnx, nny, nnz) =
     *           parts(cLcompNo)%edepo(nnx, nny, nnz) + dE
            l = l + dl
         enddo   

         dl = Move%dl - l
         cl = (l+dl/2) 
         dE = Move%dEeff*dl/Move%dl

         x = aTrack%pos%x + cl*aTrack%w%x 
         nnx = epCellIdx(x, dx)

         y = aTrack%pos%y + cl*aTrack%w%y 
         nny = epCellIdx(y, dy)

         z = aTrack%pos%z + cl*aTrack%w%z
!         nnz = sign( aint( abs(z)/dz + 1), z)  
         nnz = epCellIdx(z, dz)

         parts(cLcompNo)%edepo(nnx, nny, nnz) =
     *        parts(cLcompNo)%edepo(nnx, nny, nnz) + dE

      endif
      end 
