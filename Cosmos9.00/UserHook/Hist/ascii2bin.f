      program main
      use modHistogram1
      implicit none

      type(histogram1) h1

      integer klena
      real normf
      data normf/-1.0/     ! use normalization as already done
      integer fnoI/31/, fnoO/32/     ! file number.
      character*256 filename, Ofilename
      integer kwhistReadAscii
      integer status
      integer filel
      integer count
      integer nbinhisto
!     
!     read ascii 1-D histogram file and convert is to
!     a binary hist file
!     usage:
!       compile:   make -f ascii2bin.mk
!
!
!          get filename from command line argument
      count = NARGS()                                      

      if(count .ne. 2) then                                
         write(0,*)
     *    " ascii2bin: convert 1D ascii histogram into binary histo"
         write(0,*) "Usage: ./ascii2bin$ARCH  asciiInputFile  "
         write(0,*) " asciiInputFile: its extension must be .ahist"
         write(0,*) " If input is xxx.ahist, then output file name "
         write(0,*) " is made to be  *****  xxxx.hist **** "
         stop 100
      endif                                                

      call getarg(1, filename, status)                          
      filel = klena(filename)
      write(0,*) status, filename(1:filel)
      if(filel .le. 5 .or. 
     *    filename(filel-5:filel) .ne. ".ahist" ) then
         write(0,*) "input file name extension must be .ahist"
         stop  111
      endif
      open(fnoI, file=filename,
     *     iostat=status, access='sequential',
     *     form='formatted', action='read')
      if(status .ne. 0 ) then
         write(0,*) ' cannot open file ='
         write(0,*) filename
         stop 222
      endif
      ofilename = filename(1:filel-6)//".hist"
      open(fnoO, file=Ofilename,
     *     iostat=status, access='sequential',
     *     form='unformatted', action='write')
      if(status .ne. 0 ) then
         write(0,*) "could not open file"
         write(0,*) ofilename
         stop 333
      endif

      do while ( .true. )
!             alloc h1 and read
         nbinhisto=kwhistReadAscii(h1, fnoI)
         if(nbinhisto .le. 0 )  exit
!             not needed ?
         call kwhists(h1, normf)
!             write h1
         call kwhistw(h1, fnoO)
!             dealloc h1
         call kwhistd( h1)
      enddo
      close(fnoI)
      close(fnoO)
      end
