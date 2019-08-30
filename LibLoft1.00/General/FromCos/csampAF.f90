module modcsampAF0
  integer, parameter:: maxsampAF = 120

  type sampAF
     real(8),pointer:: x(:)
     real(8),pointer:: y(:)
     real(8),pointer:: yi(:)
     real(8),pointer:: coef(:,:)
     real(8),pointer:: coef2(:,:)
   !      coef for inverse function. may not be used
     real(8),pointer:: icoef(:,:)
   !      coef for integral upto x. may not be used  
     real(8),pointer:: igcoef(:,:)

     real(8):: sum    ! integral of function 
     integer::invf    ! becomes 1 if inverse function is used
     integer::integ   ! becomes 1 if integral from xmin to x is used
     integer:: n             ! number of points
  end type sampAF
end module modcsampAF0


module modcsampAF
  use modcsampAF0

  integer,save::sampInfoNoCounter=0

  type(sampAF),save::sampInfo(maxsampAF)

contains
  subroutine csampAF0(iowk, filen, sampInfoNo)
    implicit none
    integer,intent(in)::iowk  !  file logical number temporarily used
    character(len=*),intent(in)::filen  ! file name which contains (x,dn/dx)
    integer,intent(out)::sampInfoNo  ! the distributions defined so far
                                       ! this is sampInfoNo-th one
                                       ! so you can spcify the distribution
                   ! by sampInfoNo or directly manipulate by
                   ! sampInfo(sampInfoNo) 
!      type(sampAF),intent(out)::sampInfo  ! distribution info is stored

    integer::icon
    integer::i, n, nrow

    real(8):: x, y

    call copenf(iowk, filen, icon)
    if(icon /= 0 ) then
       write(0,*) ' error '
       write(0,*) 'file: ',filen
       write(0,*) ' could not be opened'
       stop
    endif
    i = 0
    call cskipComment(iowk, icon)
    if(icon /= 0 ) stop
!          first count the # of  data      
    do while ( .true. )
       read(iowk, *, end=100) x, y
       i  = i + 1
    enddo
100 continue
    nrow = i
! ////////////
!      write(0, *) ' nrow=', nrow
! ////////////
    if( sampInfoNoCounter >= maxsampAF ) then
       write(0,*) ' too many funcitons defined by tables > ', maxsampAF
       write(0,*) ' increase  maxsampAF in csampAF in Cosmos/cosmos'
       stop
    endif
    sampInfoNoCounter = sampInfoNoCounter +1
    sampInfoNo = sampInfoNoCounter  
    n = sampInfoNo
!   allocate memory
    allocate(sampInfo(n)%x(nrow))
    allocate(sampInfo(n)%y(nrow))
    allocate(sampInfo(n)%yi(nrow))
    allocate(sampInfo(n)%coef(nrow,1))
    allocate(sampInfo(n)%coef(nrow,2))
    allocate(sampInfo(n)%coef(nrow,3))
    allocate(sampInfo(n)%coef2(nrow,1))
    allocate(sampInfo(n)%coef2(nrow,2))
    allocate(sampInfo(n)%coef2(nrow,3))

    rewind(iowk)
    call cskipComment(iowk, icon)
    
    do i =1 , nrow
       read(iowk, *, end=200) sampInfo(n)%x(i), sampInfo(n)%y(i)
    enddo
200 continue
    close(iowk)
    sampInfo(n)%n = nrow

    call ksampAF0(sampInfo(n)%x, sampInfo(n)%y, sampInfo(n)%n, &
         sampInfo(n)%coef, sampInfo(n)%n, sampInfo(n)%yi,  &
         sampInfo(n)%sum,  sampInfo(n)%coef2 )
  end subroutine csampAF0

  subroutine csampAF0_b(iowk, filen, sampInfoNo)
!     same as csampAF0 but assumes only for getting function
!     value by interpolation. (say, for function with negative values)
    implicit none
    integer,intent(in)::iowk  !  file logical number temporarily used
    character(len=*),intent(in)::filen  ! file name which contains (x,dn/dx)
    integer,intent(out)::sampInfoNo  ! the distributions defined so far
                                       ! this is sampInfoNo-th one
                                       ! so you can spcify the distribution
                   ! by sampInfoNo or directly manipulate by
                   ! sampInfo(sampInfoNo) 
!      type(sampAF),intent(out)::sampInfo  ! distribution info is stored

    integer::icon
    integer::i, n, nrow

    real(8):: x, y

    call copenf(iowk, filen, icon)
    if(icon /= 0 ) then
       write(0,*) ' error '
       write(0,*) 'file: ',filen
       write(0,*) ' could not be opened'
       stop
    endif
    i = 0
    call cskipComment(iowk, icon)
    if(icon /= 0 ) stop
!          first count the # of  data      
    do while ( .true. )
       read(iowk, *, end=100) x, y
       i  = i + 1
    enddo
100 continue
    nrow = i
! ////////////
!      write(0, *) ' nrow=', nrow
! ////////////
    if( sampInfoNoCounter >= maxsampAF ) then
       write(0,*) ' too many funcitons defined by tables > ', maxsampAF
       write(0,*) ' increase  maxsampAF in csampAF in Cosmos/cosmos'
       stop
    endif
    sampInfoNoCounter = sampInfoNoCounter +1
    sampInfoNo = sampInfoNoCounter  
    n = sampInfoNo
!   allocate memory
    allocate(sampInfo(n)%x(nrow))
    allocate(sampInfo(n)%y(nrow))
    allocate(sampInfo(n)%yi(nrow))
    allocate(sampInfo(n)%coef(nrow,1))
    allocate(sampInfo(n)%coef(nrow,2))
    allocate(sampInfo(n)%coef(nrow,3))
!    allocate(sampInfo(n)%coef2(nrow,1))
!    allocate(sampInfo(n)%coef2(nrow,2))
!    allocate(sampInfo(n)%coef2(nrow,3))

    rewind(iowk)
    call cskipComment(iowk, icon)
    
    do i =1 , nrow
       read(iowk, *, end=200) sampInfo(n)%x(i), sampInfo(n)%y(i)
    enddo
200 continue
    close(iowk)
    sampInfo(n)%n = nrow

    call ksampAF0_b(sampInfo(n)%x, sampInfo(n)%y, sampInfo(n)%n, &
         sampInfo(n)%coef, sampInfo(n)%n)
  end subroutine csampAF0_b


  subroutine  csampAF0byArray(x, y, nrow, sampInfoNo)
!   This purpose is the same as csampAF0. the input is not
!    from file but from array
    implicit none
    integer,intent(in)::nrow  !  number of data in (x,y)
    real(8),intent(in)::x(nrow), y(nrow)  ! input data array
    integer,intent(out)::sampInfoNo  ! the distributions defined so far

    integer::n

    if( sampInfoNoCounter >= maxsampAF ) then
       write(0,*) ' too many funcitons defined by tables > ', maxsampAF
       write(0,*) ' increase  maxsampAF in csampAF in Cosmos/cosmos'
       stop
    endif
    sampInfoNoCounter = sampInfoNoCounter +1
    sampInfoNo = sampInfoNoCounter  
    n = sampInfoNo
!   allocate memory
    allocate(sampInfo(n)%x(nrow))
    allocate(sampInfo(n)%y(nrow))
    allocate(sampInfo(n)%yi(nrow))
    allocate(sampInfo(n)%coef(nrow,1))
    allocate(sampInfo(n)%coef(nrow,2))
    allocate(sampInfo(n)%coef(nrow,3))
    allocate(sampInfo(n)%coef2(nrow,1))
    allocate(sampInfo(n)%coef2(nrow,2))
    allocate(sampInfo(n)%coef2(nrow,3))

    sampInfo(n)%n = nrow
      ! copy data
    sampInfo(n)%x = x
    sampInfo(n)%y = y
    
    call ksampAF0(sampInfo(n)%x, sampInfo(n)%y, sampInfo(n)%n, &
         sampInfo(n)%coef, sampInfo(n)%n, sampInfo(n)%yi,  &
         sampInfo(n)%sum,  sampInfo(n)%coef2 )
  end subroutine csampAF0byArray

  subroutine  csampAF0byArray_b(x, y, nrow, sampInfoNo)
!   This purpose is the same as csampAF0byArray but 
!    not from file but from array and assumes
!    only for  getting function value by interpolation.

    implicit none
    integer,intent(in)::nrow  !  number of data in (x,y)
    real(8),intent(in)::x(nrow), y(nrow)  ! input data array
    integer,intent(out)::sampInfoNo  ! the distributions defined so far

    integer::n

    if( sampInfoNoCounter >= maxsampAF ) then
       write(0,*) ' too many funcitons defined by tables > ', maxsampAF
       write(0,*) ' increase  maxsampAF in csampAF in Cosmos/cosmos'
       stop
    endif
    sampInfoNoCounter = sampInfoNoCounter +1
    sampInfoNo = sampInfoNoCounter  
    n = sampInfoNo
!   allocate memory
    allocate(sampInfo(n)%x(nrow))
    allocate(sampInfo(n)%y(nrow))
    allocate(sampInfo(n)%yi(nrow))
    allocate(sampInfo(n)%coef(nrow,1))
    allocate(sampInfo(n)%coef(nrow,2))
    allocate(sampInfo(n)%coef(nrow,3))
!    allocate(sampInfo(n)%coef2(nrow,1))
!    allocate(sampInfo(n)%coef2(nrow,2))
!    allocate(sampInfo(n)%coef2(nrow,3))

    sampInfo(n)%n = nrow
      ! copy data
    sampInfo(n)%x = x
    sampInfo(n)%y = y
    
    call ksampAF0_b(sampInfo(n)%x, sampInfo(n)%y, sampInfo(n)%n, &
         sampInfo(n)%coef, sampInfo(n)%n)
  end subroutine csampAF0byArray_b


  subroutine csampAF(n, xs)
    implicit none
    integer,intent(in)::n  ! sampInfoNo 
      
    real(8),intent(out):: xs
    if( n > sampInfoNoCounter ) then
       write(0,*) ' requested AF no=',n, ' does not exist'
       stop
    endif
    call ksampAF(sampInfo(n)%x, sampInfo(n)%yi, sampInfo(n)%n, &
         sampInfo(n)%coef2, sampInfo(n)%n, xs)
  end subroutine csampAF
    
  subroutine csampAFIntp(n, xv, ans)
!         this is not for sampling but simply
!       get value of y at xv
    implicit none
    integer,intent(in):: n ! sampInfoNo obtained by csampAF0
    real(8),intent(in)::  xv  !  x value
    real(8),intent(out):: ans  ! y at xv

    if( n > sampInfoNoCounter ) then
       write(0,*) ' requested AF no=',n, ' does not exist'
       stop
    endif

    call  kcsplIntp(sampInfo(n)%x, sampInfo(n)%y, sampInfo(n)%n, &
         sampInfo(n)%coef, sampInfo(n)%n, xv, ans)
  end subroutine csampAFIntp

  subroutine csampAFmax(n,  xmax, fmax, xave)
      !       find max position and value of given function
      !      (approx value)
    implicit none
    integer,intent(in):: n ! sampInfoNo obtained by csampAF0
    real(8),intent(out):: xmax !  max position  in (x1,x2) ; approx value
    real(8),intent(out):: fmax !  max function value
    real(8),intent(out):: xave !  average value
      
    real(8):: x1, x2
    real(8):: x, dx, temp
    integer:: i

    if( n > sampInfoNoCounter ) then
       write(0,*) ' requested AF no=',n, ' does not exist'
       stop
    endif
      
    x1 = sampInfo(n)%x(1)
    x2 = sampInfo(n)%x(sampInfo(n)%n)
    x = x1
    dx = (x2-x1)/sampInfo(n)%n/10.
    xmax = x1
    call csampAFIntp(n, xmax, fmax)
    call csampAFIntp(n, x2, temp)
    if( fmax < temp ) then
       xmax = x2
       fmax = temp
    endif

    x = x + dx
    do while (x <  x2-dx/2 )
       call csampAFIntp(n, x, temp)
       if( fmax <  temp ) then
          xmax = x
          fmax = temp
       endif
       x = x + dx
    enddo

    xave = 0.
    x = x1 + dx/2
    do while (x < x2)
       call csampAFIntp(n, x, temp)
       xave = xave + temp*x*dx
       x = x + dx
    enddo
    xave = xave/sampInfo(n)%sum
 
  end subroutine csampAFmax
!      querry for the number of rows, x1, x2  of a given AF
  subroutine csampAFq(n, nrow, x1, x2, integ)
    implicit none
    integer,intent(in):: n ! sampInfoNo obtained by csampAF0
    integer,intent(out):: nrow ! number of rows of the AF
    real(8),intent(out) :: x1  ! min of x of AF
    real(8),intent(out) :: x2  ! max of x of AF
    real(8),intent(out) :: integ  ! integral of AF
    if( n > sampInfoNoCounter ) then
       write(0,*) ' csampAFqrow:'
       write(0,*) ' requested AF no=',n, ' does not exist'
       stop
    endif
    
    nrow = sampInfo(n)%n
    x1 = sampInfo(n)%x(1)
    x2 = sampInfo(n)%x(nrow)
    integ =sampInfo(n)%sum
  end subroutine csampAFq  

  subroutine csampAFatN(id, n, x, y, icon )
    implicit none
!       get (x,y) of the n-th row
    integer,intent(in):: id ! sampInfoNo obtained by csampAF0
    integer,intent(in):: n  ! specify  n-th row  
    real(8),intent(out):: x, y  !  value at n-th row.
    integer,intent(out):: icon ! =0 if n is in the ragne.
                               ! =1 if n is out of range. x,y undef.
    if( n <=0 .or. n >  sampInfo(id)%n ) then
       icon = 1
    else
       x = sampInfo(id)%x(n)
       y = sampInfo(id)%y(n)
       icon = 0
    endif
  end subroutine csampAFatN



  subroutine csampAFinvF(id, y, ans)
!     inverse function of the function specified by id
!     the function must be monotonic.
    implicit none
    integer,intent(in):: id ! sampInfoNo obtained by csampAF0
    real(8),intent(in):: y  ! argument
    real(8),intent(out):: ans  ! funciton value



    integer:: nrow ! number of rows of the AF

    nrow = sampInfo(id)%n

    if( sampInfo(id)%invf == 0 ) then
!        allocate mem. for inverse function
       allocate(sampInfo(id)%icoef(nrow,1))
       allocate(sampInfo(id)%icoef(nrow,2))
       allocate(sampInfo(id)%icoef(nrow,3))
       call kcsplCoef(sampInfo(id)%y, sampInfo(id)%x,  sampInfo(id)%n, &
            sampInfo(id)%icoef, sampInfo(id)%n)
       sampInfo(id)%invf = 1
    endif
    if( ( sampInfo(id)%y(nrow) -y)*(y-sampInfo(id)%y(1)) .lt. 0. ) then
       !  y is out of given range
       write(0,*) 'Warning: y for csampAFinvF is out of range'
       write(0,*) 'function id is', id, ' input y is', y
       write(0,*) 'function y range is ',  &
       sampInfo(id)%y(1),sampInfo(id)%y(nrow) 
    endif
    
    call  kcsplIntp(sampInfo(id)%y, sampInfo(id)%x, sampInfo(id)%n, &
         sampInfo(id)%icoef, sampInfo(id)%n, y, ans)

  end subroutine csampAFinvF


  subroutine csampAFinteg(id, x, ans) 
!         integral of the function specified by id from xmin to x
    implicit none
    integer,intent(in)::id  ! id # of the function
    real(8),intent(in)::x   !  upper bound of integral 
                            ! if it is > xmax; xmax is used
                            !          < xmin: 0 is returned
    real(8),intent(out)::ans   !  integral value

    
    integer::nrow

    nrow  = sampInfo(id)%n

    if( sampInfo(id)%integ == 0 ) then
!        allocate mem. for integral 
       allocate(sampInfo(id)%igcoef(nrow,1))
       allocate(sampInfo(id)%igcoef(nrow,2))
       allocate(sampInfo(id)%igcoef(nrow,3))
       call kcsplCoef(sampInfo(id)%x, sampInfo(id)%yi,  sampInfo(id)%n, &
            sampInfo(id)%igcoef, sampInfo(id)%n)
       sampInfo(id)%integ = 1
    endif

    if( x < sampInfo(id)%x(1) ) then
       ans  = 0.
    elseif( x > sampInfo(id)%x(nrow) ) then
       ans = sampInfo(id)%sum
    else
       call  kcsplIntp(sampInfo(id)%x, sampInfo(id)%yi, sampInfo(id)%n, &
         sampInfo(id)%igcoef, sampInfo(id)%n, x, ans)
       ans = ans* sampInfo(id)%sum
    endif
  end subroutine csampAFinteg


  subroutine  csampAFdealloc(n)
    implicit none
    integer,intent(in)::n  ! id # of the function
    
    deallocate(sampInfo(n)%x)
    deallocate(sampInfo(n)%y)
    deallocate(sampInfo(n)%yi)
    deallocate(sampInfo(n)%coef)
    deallocate(sampInfo(n)%coef2)

    if( sampInfo(n)%invf == 1 ) then 
       deallocate(sampInfo(n)%icoef)
       sampInfo(n)%invf = 0
    endif
    if( sampInfo(n)%integ == 1) then
       deallocate(sampInfo(n)%igcoef)
       sampInfo(n)%integ = 0
    endif
  end subroutine csampAFdealloc
  subroutine  csampAFdealloc_b(n)
    implicit none
    integer,intent(in)::n  ! id # of the function
    
    deallocate(sampInfo(n)%x)
    deallocate(sampInfo(n)%y)
    deallocate(sampInfo(n)%yi)
    deallocate(sampInfo(n)%coef)


    if( sampInfo(n)%invf == 1 ) then 
       deallocate(sampInfo(n)%icoef)
       sampInfo(n)%invf = 0
    endif
    if( sampInfo(n)%integ == 1) then
       deallocate(sampInfo(n)%igcoef)
       sampInfo(n)%integ = 0
    endif
  end subroutine csampAFdealloc_b

end module modcsampAF
