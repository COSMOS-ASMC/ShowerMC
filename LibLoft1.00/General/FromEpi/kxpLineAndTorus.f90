!  kxpLineAndTorus and kisInsideTorus

!   program main
!     implicit none
!     real(8)::q0(3)  ! the line passes this point
!     real(8)::a(3)   ! the direction cos of the line
!     real(8):: R=3   ! torus is assumed to be in the
!                        ! canonical coord. i.e, the central circle
!                        ! passing the torus body lies on the x-y 
!                        ! plane with the radius R.
!     real(8):: sr=2  !  torus body part has this radius (r<R)
!     
!     real(8):: k(4) ! a max of 4  crossing points are
!                       ! at q0 + k(:)a 
! 
!     real(8)::u, v, cosz, sinz, cosf, sinf
!     integer::nx, i,j, icon
!     call copenfw(20, "/tmp/kasahara/nxp",icon)
!     write(0,*) 'enter an integer <100'
!     read(*,*) j
!     do i = 1, j
!        call rndc(u)
!     enddo
!     
! !    q0(:) = (/0, 0, 0/)
! !    q0(:) = (/-3.0d0, 4.50d0, -1.5d0/)
!     q0(:) = (/-1.15474602028590d0, 3.19806363885045d0, 1.95956025169017d0/)
!!     a(:)  = (/-1.65682673662277d0, 4.04310477767844d0, 1.7756569402166d0/)
!      a(:) = (/-0.502080716336872d0, 0.845041138827991d0,-0.183903311473571d0/)
!
!
!     do while(.true.)
!        call rndc(cosz)
!        cosz =1.6*cosz-0.8
!        call kcossn(cosf, sinf)
!        sinz = sqrt(1.-cosz**2)
!!        a(:) = (/sinz*cosf, sinz*sinf, cosz/)
!        call kxpLineAndTorus(q0, a, R, sr, k, nx)
!        write(0,*) 'nx=', nx
!        write(0,*) 'k=', k(1:nx)
!        do i = 1, nx
!           if( k(i) >= 0.) then
!              write(*,*) q0
!              write(*,*)  q0(:) + k(i) * a(:)
!              write(*,*)
!              write(*,*)
!              exit
!           endif
!        enddo
!        if( nx == 0 ) then
!           write(20,*) q0
!           write(20,*) q0(:) +  7 * a(:)
!           write(20,*)
!           write(20,*)
!        endif
!        read(*,*)
!     enddo
!   end program main

    
  subroutine kxpLineAndTorus(q0, a, R, sr, k, nx)
!  get xssing point of a line and a torus; based  on
!
!　Analytic　Solution　of　intersection　between
!　　　　　　　　　High－erder　Surfaces
!
!－The　case　of　Quadric，　Fourth－order　Surfaces　（Torus）一
!
!  Yukinori　Kakazu，　Norio　OKiNo　and　Hirokazu　Watabe．
!    Bulletin　of　the　Faculty　of　Engineering，
!  　Kokkaido　University，　No．　120　（1984）

    implicit none
    real(8),intent(in)::q0(3)  ! the line passes this point
    real(8),intent(in)::a(3)   ! the direction cos of the line
    real(8),intent(in):: R     ! torus is assumed to be in the
                       ! canonical coord. i.e, the central circle
                       ! passing the torus body lies on the x-y 
                       ! plane with the radius R.
    real(8),intent(in):: sr  !  torus body part has this radius (r<R)
    
    real(8),intent(out):: k(4) ! a max of 4  crossing points are
                      ! at q0 + k(:)a 
                      ! k(1)<=k(2)...<=k(nx). 
                      ! nearest crossing point shouldb be
                      ! smallest k(i) >=0.
    integer,intent(out):: nx   !  # of real x-points. 0~4

    real(8):: b1, b2, b3, b4
    real(8):: c(0:4)
    real(8):: temp, R2
    complex(kind(0d0)):: kc(4)
    integer:: nr, ns, i

    b1 = dot_product(a, q0)*2
    b2 = dot_product(q0,q0)
    b3 = a(3)
    b4 = q0(3)

    R2 = R*R
    temp = R2 + sr**2
    c(4) = 1.0d0
    c(3) = 2*b1
    c(2) = b1**2 + 2*b2 - 2*temp + 4*R2*b3**2
    c(1) =  2*b1*b2 -2*b1*temp + 8*R2 *b3*b4
    c(0) = b2**2 -2*temp *b2 + 4*R2*b4**2 + (R2 -sr**2)**2
!       solve k^4 +C3k^3 + c2k^2 +c1k + c0 =0
    call kquarticEq(c, kc, nr, ns) 
    nx = nr
    do i = 1, nr
       k(i) = kc(i)
    enddo
  end subroutine kxpLineAndTorus

!    program main  ! test kisInsideTorus
!      implicit none
!      real(8)::R=3.
!      real(8)::sr = 2.
!      integer icon
!      real(8)::x(3)
!      integer i
!      real(8)::cosz, sinz, cosf, sinf, u
! !     call copenfw(20, "/tmp/kasahara/out",icon)
! !     do i = 1, 10000
!      do i = 1, 1
!         call rndc(cosz)
!         cosz = 2*cosz - 1.
!         sinz = sqrt(1.-cosz**2)
!         call kcossn(cosf,sinf)
!         call rndc(u)
!         u = u*5.
!         x(:) = (/sinz*cosf, sinz*sinf, cosz/)*u
! !//////
!         x(:) = (/-1.15474602028590d0, 3.19806363885045d0, &
!            1.95956025169017d0/)
! !////////////
!         call  kisInsideTorus(x, R, sr, icon)
! !///////////
!         write(0,*) 'x=',x
!         write(0,*) 'R sr =', R, sr, icon
! !///////////////
! 
!         if(icon <= 1) then
!            write(*,*)  x(:)
!         else
! !           write(20,*) x(:)
!         endif
!      enddo
!   end program main

  subroutine kisInsideTorus(x, R, sr, icon)
!      see if a given point is inside a given torus.
    implicit none
    real(8),intent(in)::x(3)  ! given point
    real(8),intent(in)::R     ! Torus circle radius. Canonical form
    real(8),intent(in)::sr    ! Torus body radius
    integer,intent(out)::icon ! 0--> inside 
                              ! 1--> (on the surface)
                              ! 2--> outside

    real(8)::R2,sr2
    real(8)::d2,f, temp

    if( abs(x(3)) > sr ) then
       icon = 2
    else
       R2 = R**2
       temp = x(1)**2 + x(2)**2
       if(  temp > (R+r)**2 ) then
          icon = 2
       elseif( temp < (R-r)**2 ) then
          icon = 2
       else
          sr2 = sr**2
          temp = sqrt(temp)
          f = (x(1)*(1. - R/temp))**2 + (x(2)*(1.-R/temp))**2 + x(3)**2 - sr2
          if( f < 0. ) then
             icon = 0
          elseif( f> 0.) then
             icon = 2
          else
             icon = 1
          endif
       endif
    end if
  end subroutine kisInsideTorus
