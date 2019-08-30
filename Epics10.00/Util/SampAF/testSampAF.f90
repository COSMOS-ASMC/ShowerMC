!       #include "csampAF.f90"  not needed
!      This is utility routine to test  random number generation
!      following an arbitrary function which is expressed in a form
!      of table.  The table format is
!      x   dn/dx
!      x   dn/dx
!      .... 
!        
!      Usage: complie by make
!
!      ./sampAF{$ARCH}  n  file-name-path 
!     where
!       n:  > 0 --->n times sampling will be tried and sampled value
!                   will be written on stdout
!             0 ---> (x,dn/dx of table with finer x step than input will be
!                   on stdout. dn/dx is normalized.
!           < 0 --> same as 0 but the data point is |n|
!       file-name-path:   the path to the file containing table
!             (table size must be maxDatPoint in Cosmos/cosmos/ZsampAF.h)
!       ./sampAF..  will write max pos, max vale, integral  of the function 
!        on stderr.
!
program main
      use modcsampAF

      implicit none
      integer :: No
      integer :: iowk  ! input file logical number temporarily used   

      character(len=256)  inputfile
      integer:: narg, lengInputFile
      integer:: n, i, nrow
        ! ny is normalized y
      real(8):: x, y, xmax, fmax, dx, xs, ny, xave, integ
      real(8):: x1, x2
      
      narg=iargc()
      if( narg .ne. 2 ) then
         write(0,*) "Usage  ./sampAFxx  n  file-name-path "
         write(0,*) "  n=0--> x, y, normalized_y @ 5xRow points  "
         write(0,*) "  n>0 -->samples x  n-times"
         write(0,*) "  n<0--> x, y, normalized_y @ |n| points"
         stop
      endif
#if defined (MacGfort) || (LinuxGfort)
      call getarg(1, inputFile)
      read(inputFile, *) n
      call getarg(2, inputFile)
      lengInputFile = len(inputFile)
#else      
      call getarg(1, inputFile, lengInputFile)
      read(inputFile, *) n
!      write(0,*) lengInputFile, inputFile
!      write(0,*) 'n=',n
      call getarg(2, inputFile, lengInputFile)
#endif      
      write(0,*) ' input file name=', inputFile(1:lengInputFile)

      iowk = 11
      call csampAF0(iowk, inputfile, No)
      call csampAFq(No, nrow, x1, x2, integ )
      write(0,*) ' # of points=', nrow
      if(n .le. 0 ) then

         if(n .lt. 0) then
            dx =( x2 - x1)/abs(n)
         else
            dx =( x2 - x1 )/nrow/5.
         endif
         x=x1
         do while (x .lt. x2)
            call csampAFIntp(No, x, y)
            write(*,'(1p,3g12.4)') x, y, y/integ
            x = x + dx
         enddo
         call csampAFIntp(No, x2, y)
         write(*,'(1p,3g12.4)')  x, y, y/integ
      else
         do i = 1, n
            call csampAF(No, xs)
            write(*,'(1p,g12.4)')  xs
         enddo
      endif

      call csampAFmax(No, xmax, fmax, xave)

      write(0,'(a, 1p, g12.4, a, g12.4)') 'x1=',x1,  'x2=', x2
      write(0,'(a, 1p, g12.4, a, g13.4,a,g12.4,a, g12.4,a,g12.4)') &
        ' max pos=', xmax, ' max valu=',fmax, ' <x>=', xave, &
        ' integ=', integ
end program 
