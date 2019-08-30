! this is to test inverse function of a function  given by table
!  and to get inegral of that funciton from xmin to x
!Usasge:
!     make  -f sampaf2.mk
!     ./a.out   n  file-name-path 
!     where
!         1 ---> x  y   table @ denser point than input (by invrse func)
!         2 --->  integral of function from xmin to x with x
!                 denser than input
!
!       file-name-path:   the path to the file containing table
!             (table size must be maxDatPoint in Cosmos/cosmos/ZsampAF.h)
!
program main
      use modcsampAF

      implicit none
      integer :: No
      integer :: iowk  ! input file logical number temporarily used   

      character(len=256)  inputfile
      integer:: narg, lengInputFile
      integer:: n, i, nrow
      real(8):: x, y, xmax, fmax, dx,  xave, integ
      real(8):: x1, x2, y1, y2, temp, ans, dy 
      
      narg=iargc()
      if( narg .ne. 2 ) then
         write(0,*) "Usage  ./a.out  n  file-name-path "
         write(0,*) "  n=1--> x y  @ 5xRow points; y is arg. x is ans "
         write(0,*) "  n=2--> interal of y form xmin to x @ 5xRow points  "
         stop
      endif
#if defined (MacGfort) || (LinuxGfort)
      call getarg(1, inputFile)
      read(inputFile, *) n
      call getarg(2, inputFile)
      lengInputFile = len(inputFile)
#else      
      call getarg(1, inputFile, lengInputFile)
!      write(0,*) lengInputFile, inputFile
      read(inputFile, *) n
!      write(0,*) 'n=',n
      call getarg(2, inputFile, lengInputFile)
#endif      
      write(0,*) ' input file name=', inputFile(1:lengInputFile)

      iowk = 11
      call csampAF0(iowk, inputfile, No)
      call csampAFq(No, nrow, x1, x2, integ )
      write(0,*) ' # of points=', nrow

      call csampAFIntp(No, x1, y1)
      call csampAFIntp(No, x2, y2)
      if( y2 < y1 ) then
         temp = y2
         y2 = y1
         y1 = temp
      endif
      dy =( y2 - y1 )/nrow/5.
      if(n == 1 ) then
         y = y1
         do while (y .lt. y2)
            call csampAFinvF(No, y, x)
            write(*,'(1p,2g12.4)') x, y
            y = y + dy
         enddo
         call csampAFInvF(No, y2, x)
         write(*,'(1p,2g12.4)')  x, y2
      elseif( n == 2 ) then
         x = x1
         dx =( x2 - x1 )/nrow/5.
         do while (x .lt. x2)
            call csampAFinteg(No, x, ans)
            write(*,'(1p,2g12.4)') x, ans
            x = x + dx
         enddo
         call csampAFinteg(No, x2, ans)
         write(*,'(1p,2g12.4)') x, ans
      else
         write(0,*) ' n =', n, ' not used' 
      endif
      call csampAFmax(No, xmax, fmax, xave)

      write(0,'(a, 1p, g12.4, a, g12.4)') 'x1=',x1,  'x2=', x2
      write(0,'(a, 1p, g12.4, a, g13.4,a,g12.4,a, g12.4,a,g12.4)') &
        ' max pos=', xmax, ' max valu=',fmax, ' <x>=', xave, &
        ' integ=', integ
end program 

