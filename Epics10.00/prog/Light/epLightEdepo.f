!      when Light is 21 or 12, energy deposit is stored in a
!      number of cells in each "light" component.
!      This is used to store the deposit in the cell.

      subroutine epLightEdepo(info, aTrack)
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



      integer::nnx, nny, nnz 
      real(8)::l, cl, hl, dl, dE      
      real(8), pointer::dx, dy, dz
      integer, pointer::nx, ny, nz
      real(8):: x, y, z
!///////////
!      call Lcompchk(" Aa ", cLcompNo)  
!////////
      dy => Lcomp(cLcompNo)%dy
      dx => Lcomp(cLcompNo)%dx
      dz => Lcomp(cLcompNo)%dz
      nx => Lcomp(cLcompNo)%nx
      ny => Lcomp(cLcompNo)%ny
      nz => Lcomp(cLcompNo)%nz

      if( .not.  allocated( parts(cLcompNo)%edepo )) then
         if( Det%cmp(Cn)%struc(1:3) == "box" ) then
            allocate( parts(cLcompNo)%edepo(1:nx,1:ny,1:nz) ) 
         if( Det%cmp(Cn)%struc(1:7) == "octagon" ) then
            allocate( parts(cLcompNo)%edepo(1:nx,1:ny,1:nz) ) 
 !             index 0 will not be used for cyl /ecyl 
         elseif( Det%cmp(Cn)%struc(1:3) == "cyl" ) then
            allocate( parts(cLcompNo)%edepo(-nx:nx, -ny:ny, 1:nz))
         elseif( Det%cmp(Cn)%struc(1:4) == "ecyl" ) then
            allocate( parts(cLcompNo)%edepo(-nx:nx, -ny:ny, 1:nz))
         elseif( Det%cmp(Cn)%struc(1:4) == "pipe" ) then
            allocate( parts(cLcompNo)%edepo(-nx:nx, -ny:ny, 1:nz))
         else
            write(0,*) ' struc=',Det%cmp(Cn)%struc(1:4) 
            write(0,*) 'not supported for epLightEdepo'
            stop
         endif
         parts(cLcompNo)%edepo = 0.
      endif
      
      if( info  ==  1 ) then
!          very short point like  track 
         nnx = aTrack%pos%x/dx + 1
         nny = aTrack%pos%y/dy + 1
         nnz = aTrack%pos%z/dz + 1
!         if(nnx .ge. 1 .and. nnx .le. nx ) then
!            if(nny .ge. 1 .and. nny .le. ny ) then
!               if(nnz .ge. 1 .and. nnz .le. nz ) then
!///////////
!      call Lcompchk(" Ab ", cLcompNo)  
!////////

         parts(cLcompNo)%edepo(nnx, nny, nnz) =
     *        parts(cLcompNo)%edepo(nnx, nny, nnz) + Move%dEeff*1.e3
!               endif
!            endif
!         endif
      else
         l = 0.
!///////////
!      call Lcompchk(" Ab ", cLcompNo)  
!////////

         dl = Lcomp(cLcompNo)%dmin 
         hl = dl/2
         dE = Move%dEeff * dl / Move%dl
         do while ( l + dl  < Move%dl )
            cl = l + hl
            x = aTrack%pos%x + cl*aTrack%w%x 
            nnx = sign( aint( abs(x)/dx + 1), x)  
            y = aTrack%pos%y + cl*aTrack%w%y 
            nny = sign( aint( abs(y)/dy + 1), y)  
            z = aTrack%pos%z + cl*aTrack%w%z
            nnz = sign( aint( abs(z)/dz + 1), z)  

            parts(cLcompNo)%edepo(nnx, nny, nnz) =
     *           parts(cLcompNo)%edepo(nnx, nny, nnz) + dE*1.e3
            l = l + dl
         enddo   

         dl = Move%dl - l
         cl = (l+dl/2) 
         dE = Move%dEeff*dl/Move%dl

         x = aTrack%pos%x + cl*aTrack%w%x 
         nnx = sign( aint( abs(x)/dx + 1), x)  

         y = aTrack%pos%y + cl*aTrack%w%y 
         nny = sign( aint( abs(y)/dy + 1), y)  

         z = aTrack%pos%z + cl*aTrack%w%z
         nnz = sign( aint( abs(z)/dz + 1), z)  

         parts(cLcompNo)%edepo(nnx, nny, nnz) =
     *        parts(cLcompNo)%edepo(nnx, nny, nnz) + dE*1.e3

      endif
      end 
