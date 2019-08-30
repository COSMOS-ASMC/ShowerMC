!
!     read binary  or ascii  lateral histogram file and
!     analyse it
!     usage:
!       compile:   make -f procLat.mk
!       see message by procLat

      implicit none
      integer binorascii, reduced, nfraca
      integer count, status
      character*80  buf
!          get filename from command line argument
      count = NARGS()                                      

      if(count .ne. 3) then                                
         write(0,*)  " ProcLat: ./procLat$ARCH",
     *   " bin   input-r.hist"
         write(0,*) " bin=1-->input is ascii file"
         write(0,*) " bin=2-->         binary file"
         write(0,*)
     *    " filename: path to the input -r.hist file(<256 chars)"
         stop 132                                        
      endif                                                
      call getarg(1, buf, status)
      read(buf,*) binorascii
      if(binorascii .ne. 1 .and. binorascii.ne. 2) then
         write(0,*) ' error input to binOrascii=',binorascii
         stop  444
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
      integer fnoT/31/     ! file number.
      character*256 filename

      integer ir, ifai
      integer code
      integer nbinhisto,  status
      integer kwhistReadAscii
      integer i
      integer  klena
      external klena
      save

      call getarg(2, filename, status)                          
!      write(0,*) ' status=',status
!      write(0,*) 'file=',filename(1:klena(filename))

      open(fnoT, file=filename,
     *         iostat=status, access='sequential',
     *         form='unformatted', action='read')
      if(status .ne. 0 ) then
         write(0,*) ' cannot open file ='
         write(0,*) filename
         stop 123
      else
!         write(0,*) filename(1:klena(filename)), 
!     *    ' opened'
      endif
      
!          web sector region
      do i = 1, ansites
         do code = 1, 4
            do ifai= 1, nfai
               read( fnoT, end=100 ) histid0
               if( histid0 .ne. '#hist1' ) then
                  write(0,*) ' histogram is not 1-D: ',histid0
                  stop 111
               endif
               call kwhistr(h10, fnoT, icon)
!                  statistical calculation
               call kwhists(h10, normf)
!                     lat analysis
               call procLat(h10,  ifai, code, i)
!                        deallocate
               call kwhistd( h10 )
            enddo
         enddo  
      enddo   
      close(fnoT)
      return
 100  continue
      write(0,*) ' unexpected EOF ' 
      stop 2345
!     **********************
      entry  procLatAscii
!     **********************
      call getarg(2, filename, status)                          
!      write(0,*) ' status=',status
!      write(0,*) 'file=',filename(1:klena(filename))

      open(fnoT, file=filename,
     *     iostat=status, access='sequential',
     *     form='formatted', action='read')
      if(status .ne. 0 ) then
         write(0,*) ' cannot open file ='
         write(0,*) filename
         stop 222
      else
!         write(0,*) filename(1:klena(filename)),
!     *      ' opened '
      endif


      do i = 1, ansites
         do code = 1, 4
            do ifai= 1, nfai
               nbinhisto=kwhistReadAscii(h10, fnoT)
               if(nbinhisto .le. 0 ) then
                  write(0,*) ' ascii read  failed'
                  stop  333
               else
!                  write(0,*) ' hist data =',nbinhisto
               endif
!               call kwhistr(h10, fnoT, icon)
!                     statistical calculation
               call kwhists(h10, normf)
!                     lat  analysis
               call procLat(h10,   ifai, code, i)
!                        deallocate
               call kwhistd( h10 )
            enddo
         enddo  
      enddo   
      close(fnoT)
      end
