!
!   Assume only 1 layer is contained in the -t.hist file
!   If multi layers are contained, need a lot of modification
!
!     read binary  or ascii  time histogram file and
!     get basic information: evid and write it in a
!     given file.
!     usage:  make -f getBasicHistoInfo.mk
!     ./getBasicHistoInfo binOrascii webonly  filename tempfile       
!
!       where binOrascii is 1 for ascii input file
!                           2 for binary //
!             webonly:  histo is from mkLDD or from FDD
!             filenname  is path to the  input histogram file
!                       (<256 char) 
!             tempfle: for memo to be used later.
      implicit none
      integer binorascii
      integer count, status
      character*80  buf

      count = NARGS()                                      
      if(count .ne. 5) then                                
         write(0,*)
     *  " must give  bin webonly filename tempfile as arguments"
         write(0,*) " bin=1-->input is ascii file"
         write(0,*) " bin=2-->         binary file"
         write(0,*)
     *   "webonly=0-->input is made from mkLDD"
     *   "  (core region data exists) "
         write(0,*)
     *   "webonly=1-->input is made from FDD rawdata"
     *   "(no core region data) "
         write(0,*) " filename: path to the input file"
         write(0,*) " tempfile: path to the temporary working file"
         write(0,*) " "
         write(0,*) " you gave ", count-1, " arguments"
         call getarg(1, buf, status)
         write(0,*) " The first one was ", buf
         stop 111                                              
      endif                                                
      call getarg(1, buf, status)
      read(buf,*) binorascii
      if(binorascii .ne. 1 .and. binorascii.ne. 2) then
         write(0,*) ' error input to binOrascii=',binorascii
         stop 111
      endif


      if(binorascii .eq. 1) then
!             ascii file
         call procTimeAscii
      else
!             bin file
         call procTimeBin
      endif
      end

      subroutine procTimeBin
      implicit none
#include "ZprivateSub.h"
#include "../../Hist/Z90histc.h"
#include "../../Hist/Z90histo.h"
#include "../../Hist/Z90hist1.h"
      character*80 buf

      type(histogram1) h10 ! 1D histogram area

      character*6 histid0  ! get histogram id  here
      integer icon
      real normf
      data normf/-1.0/     ! use normalization as already done
      integer ansites
      data ansites/1/   !  number of layers where histogram was taken
      integer fnoT/31/     ! file number for hist
      integer fnoW/32/     ! file number for working
      character*256 filename, tempfile

      integer ir, ifai
      integer code
      integer nbinhisto, mode, status, webonly
      integer kwhistReadAscii
      integer i
  

      call getarg(2, buf, status)
      read(buf, *) webonly
      call getarg(3, filename, status)
!      write(0,*) status, filename
      call getarg(4, tempfile, status)

      open(fnoT, file=filename,
     *         iostat=status, access='sequential',
     *         form='unformatted', action='read')
      if(status .ne. 0 ) then
         write(0,*) ' cannot open file ='
         write(0,*) filename
         stop 111
      endif
      open(fnoW, file=tempfile,
     *         iostat=status, access='sequential',
     *         form='formatted', action='write')
      if(status .ne. 0 ) then
         write(0,*) ' cannot open file ='
         write(0,*) tempfile
         stop 222
      endif
      
!        skip time data for core region
      if(webonly .eq. 0) then
         do i = 1, ansites
            do code = 1, 4
               read( fnoT, end=100 ) histid0
               if( histid0 .ne. '#hist1' ) then
                  write(0,*) ' histogram is not 1D: ',histid0
                  stop 111
               endif
!/////////
               write(0,*) ' reading core region'
!/////////////
               call kwhistr(h10, fnoT, icon)
               call kwhists(h10, normf)
!                   deallocate
               call kwhistd( h10 )
            enddo
         enddo
      endif
!          web sector region
!           open memo file
!//////////////
      write(0,*) ' reading web region'
!//////////////
      do i = 1, ansites
         do code = 1, 4
            do ifai= 1, nfai
               do ir= 1, nrbin
                  read( fnoT, end=100 ) histid0
                  if( histid0 .ne. '#hist1' ) then
                     write(0,*) ' histogram is not 1-D: ',histid0
                     stop 111
                  endif
                  call kwhistr(h10, fnoT, icon)
                  write(fnoW,'(a)') h10%c%id
!                        deallocate
                  call kwhistd( h10 )
                  close(fnoW)
                  close(fnoT)
                  return  !  *********
               enddo
            enddo  
         enddo   
      enddo
      return
 100  continue
      write(0,*) ' unexpected EOF ' 
      stop 2345
!     **********************
      entry  procTimeAscii
!     **********************

      call getarg(2, buf, status)
      read(buf, *)  webonly
      call getarg(3, filename, status)                          
!      write(0,*) status, filename
      call getarg(4, tempfile, status)                          
      open(fnoT, file=filename,
     *     iostat=status, access='sequential',
     *     form='formatted', action='read')
      if(status .ne. 0 ) then
         write(0,*) ' cannot open file ='
         write(0,*) filename
         stop 111
      endif
      open(fnoW, file=tempfile,
     *     iostat=status, access='sequential',
     *     form='formatted', action='write')
      if(status .ne. 0 ) then
         write(0,*) ' cannot open file ='
         write(0,*) tempfile
         stop 111
      endif

      if(webonly .eq. 0 ) then
!        skip time data for core region
         do i = 1, ansites
            do code = 1, 4
               nbinhisto=kwhistReadAscii(h10, fnoT)
               if(nbinhisto .le. 0 ) then
                  write(0,*) ' ascii read  failed'
                  stop 111
               endif
               call kwhists(h10, normf)

!     deallocate
               call kwhistd( h10 )
            enddo
         enddo
      endif
!          web sector region
      do i = 1, ansites
         do code = 1, 4
            do ifai= 1, nfai
               do ir= 1, nrbin
                  nbinhisto=kwhistReadAscii(h10, fnoT)
                  if(nbinhisto .le. 0 ) then
                     write(0,*) ' ascii read  failed'
                     stop 111
                  endif
                  write(fnoW,'(a)') h10%c%id
!                        deallocate
                  call kwhistd( h10 )
                  close(fnoW)
                  close(fnoT)
                  return  !  *********
               enddo
            enddo  
         enddo   
      enddo
      end
