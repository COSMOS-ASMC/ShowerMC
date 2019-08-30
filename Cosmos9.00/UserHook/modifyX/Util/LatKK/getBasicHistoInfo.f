!  This is diff. from the one for "time".
!   Assume only 1 layer is contained in the -r.hist file
!   If multi layers are contained, need a lot of modification
!
!     read binary  or ascii  r histogram file and
!     get basic information: evid and write it in a
!     given file.
!     usage:  make -f getBasicHistoInfo.mk
!     ./getBasicHistoInfo binOrascii  filename tempfile       
!
!       where binOrascii is 1 for ascii input file
!                           2 for binary //
!             filenname  is path to the  input histogram file
!                       (<256 char) 
!             tempfle: for memo to be used later.
      implicit none
      integer binorascii
      integer count, status
      character*80  buf

      count = NARGS()                                      
      if(count .ne. 4) then                                
         write(0,*)
     *  " must give  bin filename tempfile as arguments"
         write(0,*) " bin=1-->input is ascii file"
         write(0,*) " bin=2-->         binary file"
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
         call procLatAscii
      else
!             bin file
         call procLatBin
      endif
      end

      subroutine procLatBin
      implicit none
#include "../ZprivateSub.h"
#include "../../../Hist/Z90histc.h"
#include "../../../Hist/Z90histo.h"
#include "../../../Hist/Z90hist1.h"
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
      integer nbinhisto, mode, status
      integer kwhistReadAscii
      integer i
  

      call getarg(2, filename, status)
!      write(0,*) status, filename
      call getarg(3, tempfile, status)

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

      do i = 1, ansites
         do code = 1, 4
            do ifai= 1, nfai
               read( fnoT, end=100 ) histid0
               if( histid0 .ne. '#hist1' ) then
                  write(0,*) ' histogram is not 1-D: ',histid0
                  stop 111
               endif
               call kwhistr(h10, fnoT, icon)
               write(fnoW,'(a)') h10%c%id
!     deallocate
               call kwhistd( h10 )
               close(fnoW)
               close(fnoT)
               return           !  *********
            enddo
         enddo  
      enddo
      return
 100  continue
      write(0,*) ' unexpected EOF ' 
      stop 2345
!     **********************
      entry  procLatAscii
!     **********************

      call getarg(2, filename, status)                          
!      write(0,*) status, filename
      call getarg(3, tempfile, status)                          
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

!          web sector region
      do i = 1, ansites
         do code = 1, 4
            do ifai= 1, nfai
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
               return           !  *********
            enddo
         enddo  
      enddo   
      end
