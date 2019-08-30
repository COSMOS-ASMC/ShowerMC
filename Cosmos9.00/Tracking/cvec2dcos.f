!       cvec2dcos      :  2D vector angle to 3D direction cos 
!       cdispVec2cos   :  2D displacement vector to 3D directon cos and length
!
      subroutine cvec2dcos(vx, vy, dircos)
!        small vector( vx, vy ) is converted into direction cos.
!        valid only for samll vec, which is alw.r(2)s assumed.
      implicit none

#include  "Zcoord.h"
      real*8 vx, vy
      type(coord)::dircos
!
      real*8 sum
!
      if(vx .eq. 0 .and. vy .eq. 0.) then
         dircos%r(1) = 0. 
         dircos%r(2) = 0.
         dircos%r(3) = 1.
      else
         if(abs(vx) .lt. 1.) then
            dircos%r(1) = vx
         else
            dircos%r(1) = 1.
         endif
         if(abs(vy) .lt. 1.) then
            dircos%r(2) = vy
         else
            dircos%r(2) = 1.
         endif
         dircos%r(3) = max(0.d0, 1.d0 - (vx**2 + vy**2)/2)
!               make consistent as direction cos
         sum = sqrt(dircos%r(1)**2 + dircos%r(2)**2 + dircos%r(3)**2)
         dircos%r(1) = dircos%r(1)/sum
         dircos%r(2) = dircos%r(2)/sum
         dircos%r(3) = dircos%r(3)/sum
      endif
      end
!          displacement vector is converted into direction cos
!          len is the length of vector.
      subroutine cdispVec2Dcos(dx, dy, len, dircos)
!          len:  output
!        dircos: output
!
      implicit none
!----      include 'Zcoord.h'
#include  "Zcoord.h"
      real*8 dx, dy, len
      type(coord)::dircos
!
      if(dx .eq. 0. .and.  dy .eq. 0.) then
         len = 0.
         dircos%r(1) = 0.
         dircos%r(2) = 0
         dircos%r(3) = 1.
      else
         len = sqrt(dx**2 + dy**2)
         dircos%r(1) = dx/len
         dircos%r(2) = dy/len
         dircos%r(3) = 0.
      endif
      end

       


         
         
      
