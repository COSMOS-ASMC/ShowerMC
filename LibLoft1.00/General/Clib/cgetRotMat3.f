!c       to test cgetRotMat3, cinvRotMat3, cmultRotMat3, capplyRot3
!c           cgetRotMat3: get rotation matrix of x/y/z axis
!c           cinvRotMat3: invert a rotation matrix
!c           cmultRotMat3:  multiply two rotation matrix
!c           capplyRot3: applay a rotation on a 3d vector.
!c
!         implicit none
!         real*8 rm(3, 3),  rmn(3,3), rmm(3, 3)
!         real*8  Torad, vn(3), v(3)
!         integer i, j, m
!         v(1) =1.
!         v(2) =0.
!         v(3) = 1.
!
!         do i=1, 3
!            do j=1, 3
!                rm(i, j)=1.d30
!                rmn(i, j)=-1.d30
!                rmm(i, j)=1.d15
!            enddo
!         enddo
!         do i=1, 3
!            vn(i)= 1.e50
!         enddo
!         Torad=asin(1.d0)/90.d0
!         m=3
!         call cgetRotMat3(m, 1.d0,  0.d0, rm)
!         write(*,*) rm
!         call cgetRotMat3(m, cos(45.*Torad), -sin(45.*Torad), rm)
!         call cinvRotMat3(rm, rmn)
!         call cgetRotMat3(m,cos(-45.*Torad), -sin(-45.*Torad), rmm)
!         write(*, *) ' rotation matrix around m=', m
!         do i=1, 3
!             write(*,*) (rm(i,j),j=1, 3)
!          enddo
!         write(*,*) "-----inverse of the above--"
!         do i=1, 3
!             write(*,*) (rmn(i,j),j=1, 3)
!          enddo
!         write(*,*) "-----should be the same as above--"
!         do i=1, 3
!             write(*,*) (rmm(i,j),j=1, 3)
!          enddo
!         write(*,*) "-------"
!         call capplyRot3(rm, v, vn)
!         write(*, *) ' v is rotated by rm'
!         write(*,*) vn
!         call cmultRotMat3(rm, rmn, rmm)
!         write(*, *)' unit matrix should be seen'
!         write(*, *) rmm
!         call capplyRot3(rmm, v, vn)
!         write(*, *) " no change from v"
!         write(*,*) vn
!        end
!        **************************************************************
!        *
!        *  cgetRotMat3: make rotation matrix of 3-d coordinate axes.
!        *
!        **************************************************************
!  /usage/  call cgetRotMat3(m, cosa, sina, rm)
!
!    m: input. integer  1--> rotaion of coordinate around x axis
!                       2-->   //                         y
!                       3-->   //                         z
!                     rotation is made anticlock wise (e.g,
!                     if m=3,
!                                  !y    / new.r(1)
!                    $  new.r(2)      !    /
!                        $         !   / *
!                            $     !  /   * ang ( if > 0)
!                                $ ! /     *
!                                   ~~~~~~~~~~~~~~~~~ x
!                         z is directed to the eyes.
!   cosa: input. real*8. cos of  rotation angle.
!   sina: input. //      sin //
!   rm(3,3): output. real*8  rotation matrix.  v'=rm*v is the new coordinate
!                   of a vector v(see capplyRot3)
!
      subroutine cgetRotMat3(m, cosa, sina, rm)
          implicit none
!
          integer m
          real*8 cosa, sina
          real*8 rm(3, 3)
!
          integer i, j, m1, m2
          character*70 msg
!
          if(m .ge. 1 .and. m .le. 3) then
             do i=1, 3
                do j=1, 3
                   rm(i, j)= 0.d0
                enddo
             enddo   
             rm(m, m)=1.d0
             m1=mod(m,3)+1
             m2=mod(m1,3)+1
             rm(m1,m1)=cosa
             rm(m2,m2)=cosa
             rm(m1,m2)=sina
             rm(m2,m1)=-sina
          else
             write(msg, *) ' invalid m=',m,' to cgetRotMat3 '
             call cerrorMsg(msg,0)
          endif
      end
      subroutine cinvRotMat3(rm, rmn)
!             Invert rotation matrix rm and put into rmn.
!             rm should be a roation matrix made by calling
!             cgetRotMat3 (with cosa, sina).  rmn can be made by calling
!             cgetRotMat3 with cosa, -sina, too.  
!             rmn cannot be the same arrays as rm.
!             rmn is nothing but the transposed matrix of rm.
          implicit none
          real*8 rm(3,3), rmn(3,3)
          integer i, j
          do  i=1,3
             do  j=i+1, 3
                rmn(i,j)=rm(j, i)
                rmn(j,i)=rm(i, j)
             enddo
          enddo
          do  i=1, 3
             rmn(i,i)=rm(i, i)
          enddo
       end
       subroutine cmultRotMat3(a, b, c)
!            3-d matrix product c=a*b
!            c cannot be either of a or b.
          implicit none
          real*8 a(3,3), b(3,3), c(3,3)
          integer i, j, k
          real*8 ab
          do  i=1, 3
             do j=1, 3
                ab=0.
                do k=1, 3
                   ab=ab+ a(i,k)*b(k,j)
                enddo
                c(i,j)=ab
             enddo   
          enddo
       end
       subroutine capplyRot3(a, v, vn)
!          3-d transformation matrix a is multiplied by
!          a vector v to obtain a new vector vn.
!          vn can  be v.
!
           implicit none
           real*8 a(3,3), v(3), vn(3)
           
!
           real*8 sum
           integer i, j
!
          do  i=1, 3
             sum=0.
             do  j=1, 3
                 sum=sum + a(i, j)*v(j)
             enddo    
             vn(i)=sum
          enddo
       end
