      subroutine xsecfnc(npar,gin,f,paramx,iflag)
      use modxsecFit 
      implicit none
! npar: input. number of current variable  parameters
! gin:  optional output of gradient
!  f:  output. function value to be minimized 
! paramx: input. vector of const and variable parameters
!  iflag: input. depending on this value, what to do 
!        is determined.
!       1--> some preparation for computing f
!       2--> compute gin
!       3--> fittng finished
!       for all other cases: we  must compute f
!
      integer  npar
      real*8   gin(*)
      real*8   f
      real*8   paramx(*)
      integer iflag
!
      integer i
      real*8 fval,  xx, temp
      real*8  a, b, c
      save
      if (iflag .eq. 1) then
!         write(0,*) 'current npar=',npar,  ' npoint=',npoint
!         write(0,*) ' param=',(paramx(i), i=1,npar)
      endif
      if(iflag .eq. 2 ) then
         write(0,*) ' no grad computed'
         stop
      endif
!        compute f

      chisq = 0.
      do  i= 1, npoint
         a = paramx(1)
         b =  paramx(3)
         c = paramx(4)
         temp = log(x(i))
         fval =  a + (b+c*temp)* temp
         chisq= chisq + (fval-y(i))**2/y(i)
      enddo
      f = chisq

      if(iflag .eq. 3) then
!         do i = 1, npoint
!            fval =   paramx(1)*x(i)**(paramx(2) +
!            write(*,*) x(i), y(i), fval
!         enddo
         do i = 1, npar
            oparam(i) = paramx(i)
         enddo
      endif
      end
