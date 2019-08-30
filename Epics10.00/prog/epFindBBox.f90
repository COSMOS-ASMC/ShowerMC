 module FindBBox
   implicit none
   real(8),save:: Ix, Iy, Iz, Ixy, Iyz, Izx
   real(8),save:: I(3,3)
   real(8),save:: abc(3)
   real(8),save:: Eigen(3)
   real(8),save:: EigenVec(3,3), conv(3,3)
   real(8),save:: minmax(3,2)
   real(8),save:: rg(3)
   real(8),save:: dig=10000.0d0

   contains

 subroutine epFindBBox(rin, n,  abcout, rm,  org)
!   Find minimal  bounding box which contains all the given points.
!   Minimual here is not mathematially proven but fairly good in 
!   many practical cases. (A mathematical rigorous treatment seems
!   to exist--Joseph O'Rourke: 1985), but it is complex and takes 
!   O(n^3) time, and it is said some heuristic approaches are still
!   being tried. 
!     The method here is hinted by the moment of inertia  matrix
!   which gives 3 principal axes. The reference frame could 
!   be formed by using such axes as x-,y-,z- axis. (We call canonical
!   system, or frame or space).  The min and max of x,y,z there will be
!   used to form a box.
!     However, if we treat each input point as having mass 1, the principal
!   axes sometimes go to undesired direction. The moment of inertia matrix
!   uses square of length, so we use length propotional quantity instead.
!   For example, Ix= sum (y^2 + z^2) is used in the moment of inertia matrix
!   but we use Ix = sqrt( sum(y^2+z^2)). This change seems to work very good.
!
!   Strategy: 
!    1) Get c.o.g (=rg)(mass of each point is 1) of the given points, r.
!    2) Do conversion: r0 = r-rg
!    3) Get 'inertia matrix' I  from  r0.
!    4) Get eigenvalue and eigenvector (matrix) of I: EigenVec
!    5) Get inverse of eigenvector; conv  ( conv * EigenVec = Unit matrix)
!    6) Convert r0 by  applying conv to r0. (canonical space)
!    7) Get min and max in the cannonical space.  
!    8) Put little margin to min and max.
!    9) Return abc, rotation matrix (rm=EigenVec), and origin.
!       to be useable for definition of inclined box. 
! Auxiliary subroutine:
! After calling epFindBBox, for the current rin
!      real(8):: rg(3), conv(3,3), bboxc(3,8), bboxo(3,8)
!        call epqFindBBox(rg, conv, bboxc, bboxo)
!     rg:  c.o.g
!     conv: inverse of EigenVec.  veco(:)=matmul(conv, vecc(:)) + org(:)
!           will convert vecc in canonical space into original space
!    bboxc: bounding box vertex in canonical space. 
!    bboxo: //                     original  spaec. 
!      first 4 and last four are coplanar, respectively.
! For other variable,  we can use
! module FindBBOX to get them directly.
!

   implicit none
   integer,intent(in):: n  ! # of points
   real(8),intent(in):: rin(3,n)  ! input points (3 is x,y,z)
   real(8),intent(out):: abcout(3)   ! bounding box edge length, a, b, c 
   real(8),intent(out):: rm(3,3)  ! =EigenVec. canonical space to c.o.g space conv
                   !  is possilbe via this matrix.
                   !  rm(1:3,1) and rm(1:3,2) are  the same  as the direction
                   !  cosines of box's x  and y axis which are used in
                   !  Epics config file to denote rotation of a given box: e.g
                   !  10 box Pb  0 0 0 / ox oy oz  a b c  Xx Xy Xz Yx Yy Yz
                   !  (Xx,Xy,Xz) is rm(1:3,1).  (Yx, Yy, Yz) = rm(1:3,2).
                   !  (ox, oy, oz) is org(3) below.  rm(1:3,3) is automatically 
                   !  fixed since rm is orthogonal.
                   !  r in the canonical space can be converted to the one in
                   !  c.o.g. space by matmul(rm, r).

   real(8),intent(out):: org(3)
   real(8),allocatable:: r(:,:), r0(:,:)
   integer::j
   integer:: iflag(3)
   real(8):: work(3,6)
   integer:: iwk(3)
   integer:: ierr
   integer:: indx(3)
   real(8),save:: eps= 1.d-8
!#if defined (BBOXDEBUG)
!   real(8):: unit(3,3)
!#endif

   allocate( r(3,n) )
   allocate( r0(3,n) )

   do j = 1, 3
      rg(j) = sum( rin(j,1:n)) /n
      r(j,1:n) = rin(j,1:n) - rg(j)
   enddo

!     inertial moment is ^2 and not suitable so
!   we take sqrt instead
   Ix = sum( r(2,1:n)**2 + r(3,1:n)**2)
   Ix = sqrt(IX) !!!

   Iy = sum( r(3,1:n)**2 + r(1,1:n)**2)
   Iy = sqrt(Iy) !!!

   Iz = sum( r(1,1:n)**2 + r(2,1:n)**2)
   Iz = sqrt(Iz)  !!!

   Ixy = - sum(r(1,1:n) * r(2,1:n) ) ! put - here.
   Ixy = sign( sqrt(abs(Ixy)), Ixy)  !!!

   Iyz = - sum(r(2,1:n) * r(3,1:n) )
   Iyz = sign( sqrt(abs(Iyz)), Iyz)  !!!

   Izx = - sum(r(3,1:n) * r(1,1:n) )
   Izx = sign( sqrt(abs(Izx)), Izx)  !!!
!    matrix;  this I is same as i. so you cannot use integer i.   
   I(1,1) = Ix
   I(1,2) = Ixy
   I(2,1) = Ixy
   I(2,2) = Iy
   I(1,3) = Izx
   I(3,1) = Izx
   I(3,3) = Iz
   I(2,3) = Iyz 
   I(3,2) = Iyz 
!        next is from HitachSSL. with debug/check mode it will collapse 
!        so don't use it.
   ! call khf2m(I, 3, 3, 3, 3, eps, 2, Eigen, EigenVec, iflag, &
   !      work, iwk,  ierr)
!        next is from Maruzen SSL. need some bug fix.
   call eigrs(I, 3, 3, 3, 3, eps, work, iwk, Eigen, EigenVec, ierr)

   !          sort Eigen / EigenVec  as ascending order. needed if
   !  eigrs is used; in khf2m case, already sorted
   !   indx gets only index for ascending order
   call  kqsortd(Eigen,indx, 3)
   do j = 1, 3
      work(j, 1) = Eigen(indx(j))  ! move eigen value to work
   enddo
   do j = 1, 3
      Eigen(j) = work(j,1)   !  restore them
   enddo
   do j = 1, 3
      work(:, j) = EigenVec(:,indx(j))  ! same for eigenvector
   enddo
   do j = 1, 3
      EigenVec(:,j) = work(:,j)    ! restore
   enddo
!          end of sort

!     EigenVec is the matrix to convert (x,y,z) in canonical space
!     into c.o.g system. (The canonical system is such that I becomes
!     diagonal. 
!     conv:  matrix to convert  (x,y,z) in  c.o.g system into 
!          canonical system where EigenVec forms x,y,z axes.
!          it is the inverse of EigenVec and since I is symmetric,
!          it is just transpose of EigenVec (orthogonal)

   rm(:,:) = EigenVec(:,:)
   conv(:,:) = transpose(EigenVec(:,:))
   r0(:,1:n) = matmul(conv(:,:), r(:,1:n))   ! r0: canonical 

   do j = 1, 3
      minmax(j, 1) = minval( r0(j,:)  )
      minmax(j, 2) = maxval( r0(j,:)  )
   enddo
!      abc= max-min but  put some margin (see below)
!   dig=10000. case. 
!         x            int(dig*x)/dig   // + 1./dig      //-1./dig
!   3.14159265358979  3.14150000000000  3.14160000000000  3.14140000000000     
!  -3.14159265358979 -3.14150000000000 -3.14140000000000 -3.14160000000000     
!         0                 0          1.0000000000E-004 -1.000000000E-004
   minmax(:,1) = int(dig*minmax(:,1)-2.0d0)/dig  ! min -2 for safety of org
   minmax(:,2) = int(dig*minmax(:,2)+1.0d0)/dig  ! max
   org(:) =  matmul(EigenVec(:,:), minmax(:,1))  + rg(:)  !
   org(:) = int(dig*org(:))/dig
   abc(:) =minmax(:,2) - minmax(:,1)    ! get abc
   abcout(:) = abc(:)

!#if defined (BBOXDEBUG)
!   write(0,*) ' c.o.g'
!   write(0,'(3g15.5)') org(:)
!   write(0,*) ' ierr=',ierr
!   ! write(0,*) ' iflag =', iflag(:)
!   write(0,*) 'eigen value'
!   write(0,'(3g15.5)') Eigen(:)
!!            This type of output may result in a warning
!!            if ifort is used with -check all
!   write(0,*) ' eigen vector'
!   write(0,'(3g15.5)')  EigenVec(1,:)
!   write(0,'(3g15.5)')  EigenVec(2,:)
!   write(0,'(3g15.5)')  EigenVec(3,:)
!   write(0,*) ' eigen vector for drawing'
!   write(0,'(3g15.5)')  EigenVec(:,1)*10
!   write(0,'(3g15.5)')  0, 0, 0
!   write(0,*)
!   write(0,*)
!   write(0,'(3g15.5)')  EigenVec(:,2)*7
!   write(0,'(3g15.5)')  0, 0, 0
!   write(0,*)
!   write(0,*)
!   write(0,'(3g15.5)')  EigenVec(:,3)*5
!   write(0,'(3g15.5)')  0, 0, 0
!   write(0,*)
!   write(0,*)
!   write(0,*) ' Transposed EigenVec: conv'
!   write(0,'(3g15.5)') conv(1,:)
!   write(0,'(3g15.5)') conv(2,:)
!   write(0,'(3g15.5)') conv(3,:)
!
!   unit(:,:) = matmul(EigenVec, conv)
!   write(0,*) ' EigenVec x conv should be unit mat'
!   write(0,'(3g15.5)')  unit(1,:)
!   write(0,'(3g15.5)')  unit(2,:)
!   write(0,'(3g15.5)')  unit(3,:)
!#endif

   deallocate( r )
   deallocate( r0 )
 end subroutine epFindBBox


 subroutine epqFindBBox(rgout, convout, bboxvcout, bboxvoout)
  implicit none
  real(8),intent(out) :: rgout(3)
  real(8),intent(out) :: convout(3,3)
  real(8),intent(out) :: bboxvcout(3,8)
  real(8),intent(out) :: bboxvoout(3,8)
  
  integer::j

  rgout(:) = rg(:)
  convout(:,:) = conv(:,:)

!     in the canonical space.
  bboxvcout(:,1) = minmax(:,1)     ! left bottom vertex
  bboxvcout(3,1:4) = minmax(3, 1)  ! z values of bottom surface

  bboxvcout(1,  2) = bboxvcout(1,1) + abc(1)    
  bboxvcout(2,  2) = bboxvcout(2,1) 

  bboxvcout(1,3)  =  bboxvcout(1,2)
  bboxvcout(2,3)  =  bboxvcout(2,2) + abc(2)

  bboxvcout(1,4)  = bboxvcout(1,1)
  bboxvcout(2,4)  = bboxvcout(2,3)

  bboxvcout(:,5:8) = bboxvcout(:,1:4)
  bboxvcout(3,5:8) = bboxvcout(3,1:4) + abc(3)
  do j = 1, 8
     bboxvoout(:,j) = matmul(EigenVec(:,:), bboxvcout(:,j)) + rg(:)
  enddo
end subroutine epqFindBBox
end module FindBBox
