!c       to test cgetRotMat4, cinvRotMat4, cmultRotMat4, capplyRot4
!c           cgetRotMat4: get rotation matrix of x/y/z axis
!c           cinvRotMat4: invert a rotation matrix
!c           cmultRotMat4:  multiply two rotation matrix
!c           capplyRot4: applay a rotation on a 3d vector.
!c
!         implicit none
!         include '../Zptcl.h'
!         type(fmom):: vn, v
!         real*8 rm(4, 4),  rmn(4,4), rmm(4,4)
!         real*8  Torad
!         integer i, j, m
!         v.x =1.
!         v.y =0.
!         v.z =1.
!c
!         do i=1, 4
!            do j=1, 4
!                rm(i, j)=1.d30
!                rmn(i, j)=-1.d30
!                rmm(i, j)=1.d15
!            enddo
!         enddo
!         do i=1, 3
!            vn.p(i)= 1.e50
!         enddo
!         Torad=asin(1.d0)/90.d0
!         m=3
!         call cgetRotMat4(m, 45.*Torad, rm)
!         call cinvRotMat4(rm, rmn)
!         call cgetRotMat4(m, -45.*Torad, rmm)
!         write(*, *) ' rotation matrix around m=', m
!         do 10 i=1, 4
!             write(*,*) (rm(i,j),j=1, 4)
!  10     continue
!         write(*,*) "-----inverse of the above--"
!         do 20 i=1, 4
!             write(*,*) (rmn(i,j),j=1, 4)
!  20     continue
!         write(*,*) "-----should be the same as above--"
!         do 30 i=1, 4
!             write(*,*) (rmm(i,j),j=1, 4)
!  30     continue
!         write(*,*) "-------"
!         call capplyRot4(rm, v, vn)
!         write(*, *) ' v is rotated by rm'
!         write(*,*) vn.p
!         call cmultRotMat4(rm, rmn, rmm)
!         write(*, *)' unit matrix should be seen'
!         write(*, *) rmm
!         call capplyRot4(rmm, v, vn)
!         write(*, *) " no change from v"
!         write(*,*) vn.p
!        end
!        **************************************************************
!        *
!        *  cgetRotMat4: make rotation matrix of 3-d coordinate axes.
!        *
!        **************************************************************
!  /usage/  call cgetRotMat4(m, ang, rm)
!
!    m: input. integer  1--> rotaion of coordinate around x axis
!                       2-->   //                         y
!                       3-->   //                         z
!                     rotation is made anticlock wise (e.g,
!                     if m=3,
!                                  !y    / new x
!                    $  new y      !    /
!                        $         !   / *
!                            $     !  /   * ang ( if > 0)
!                                $ ! /     *
!                                   ~~~~~~~~~~~~~~~~~ x
!                         z is directed to the eyes.
!       ang: input. real*8.  rotation angle in radian.
!   rm(4,4): output. real*8  rotation matrix.  v'=rm*v is the new coordinate
!                   of a vector v(see capplyRot4)
!
      subroutine cgetRotMat4(m, ang, rm)
          implicit none
!
          integer m
          real*8 ang
          real*8 rm(4, 4)
!
          integer i, j, m1, m2
          real*8 c, s
!
          if(m .ge. 1 .and. m .le. 3) then
             do i=1, 4
                do j=1, 4
                   rm(i, j)= 0.d0
                enddo
             enddo   
             rm(4,4)=1.d0
             rm(m, m)=1.d0
             m1=mod(m,3)+1
             m2=mod(m1,3)+1
             c=cos(ang)
             s=sin(ang)
             rm(m1,m1)=c
             rm(m2,m2)=c
             rm(m1,m2)=s
             rm(m2,m1)=-s
          else
             write(*,*) ' invalid m=',m,' to cgetRotMat4 '
             stop  9999
          endif
      end
      subroutine cinvRotMat4(rm, rmn)
!             Invert rotation matrix rm and put into rmn.
!             rm should be a roation matrix made by calling
!             cgetRotMat4 (with ang).  rmn can be made by calling
!             cgetRotMat4 with -ang, too.  this one is to avoid
!             computing the cos and sin for reducing time.
!             rmn cannot be the same arrays as rm.
!             rmn is nothing but the transposed matrix of rm.
          implicit none
          real*8 rm(4,4), rmn(4,4)
          integer i, j
          do  i=1,4
             do  j=i+1, 4
                rmn(i,j)=rm(j, i)
                rmn(j,i)=rm(i, j)
             enddo
          enddo
          do  i=1, 4
             rmn(i,i)=rm(i, i)
          enddo
       end
       subroutine cmultRotMat4(a, b, c)
!            3-d matrix product c=a*b
!            c cannot be either of a or b.
          implicit none
          real*8 a(4,4), b(4,4), c(4,4)
          integer i, j, k
          real*8 ab
          do  i=1, 4
             do j=1, 4
                ab=0.
                do k=1, 4
                   ab=ab+ a(i,k)*b(k,j)
                enddo
                c(i,j)=ab
             enddo   
          enddo
       end
       subroutine capplyRot4(a, v, vn)
!          3-d transformation matrix a is multiplied by
!          a vector v to obtain a new vector vn.
!          vn can  be v.
!
           implicit none
!----           include '../Zptcl.h'
#include  "Zptcl.h"
           type(fmom):: v, vn
           real*8 a(4,4)
           type(fmom):: tmp
!
           real*8 sum
           integer i, j
!
          do  i=1, 3
             sum=0.
             do  j=1, 4
                 sum=sum + a(i, j)*v%p(j)
             enddo    
             tmp%p(i)=sum
          enddo
          vn = tmp
       end

