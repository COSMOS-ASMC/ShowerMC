  subroutine kSRunge_Kutta4(t, u, nf, nsim,  h, f)
    implicit none
    integer,intent(in):: nf  ! t(0), t(1),...t(nf)
                       ! contains t values with step h
    integer,intent(in):: nsim  ! nsim simultaneous equations
              
    real(8),intent(inout)::t(0:nf), u(nsim, 0:nf)
    real(8),intent(in):: h
    interface
       function f(t,u,nsim) result(array)
         implicit none
         real(8),intent(in)::t
         integer,intent(in)::nsim
         real(8),intent(in)::u(nsim)
         real(8)::array(nsim)
       end function f
    end interface
    integer::n, i
    real(8),allocatable:: u0(:),  k1(:), ks(:), k2(:), k3(:), k4(:), ux(:)

    real(8):: hh

    allocate(u0(nsim))
    allocate(k1(nsim))
    allocate(k2(nsim))
    allocate(ks(nsim))
    allocate(k3(nsim))
    allocate(k4(nsim))
    allocate(ux(nsim))
    hh = h/2
    n = 0
    do while ( n < nf)
       u0(:) = u(:,n)
       k1(:) = h*f(t(n), u0, nsim)

       ux(:) = u0(:) + k1(:)/2
       k2(:) =  h*f(t(n)+hh, ux,nsim)

       ux(:) =  u0(:)+ (k1(:)+k2(:))/4.d0
       ks(:) = h*f(t(n)+hh, ux, nsim)

       ux(:) =  u0(:)+ks(:)/2
       k3(:) = h*f(t(n)+hh, ux, nsim)


       ux(:) =   u0(:)+ks(:)
       k4(:) = h*f(t(n)+h, ux, nsim)

       u(:,n+1) = u0(:) + (k1(:) + (k2(:)+ k3(:))*2 + k4(:))/6.d0
       n = n + 1
    enddo
    deallocate(u0)
    deallocate(k1)
    deallocate(k2)
    deallocate(ks)
    deallocate(k3)
    deallocate(k4)
    deallocate(ux)
  end subroutine kSRunge_Kutta4

       
