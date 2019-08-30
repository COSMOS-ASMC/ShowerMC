!
!     read binary  or ascii  time histogram file and
!     analyse it
!     usage:
!       compile:   make -f procTime.mk
!       see message by procTime

      implicit none
#include "ZprivateSub.h"
      integer binorascii, reduced, nfraca
      integer count, status
      character*80  buf
!          get filename from command line argument
      count = NARGS()                                      

      if(count .ne. 7) then                                
         write(0,*)  " ProcTime: ./procTime$ARCH",
     *   " bin reduced nfraca  webonly smooth input-t.hist"
         write(0,*) " bin=1-->input is ascii file"
         write(0,*) " bin=2-->         binary file"
         write(0,*) " reduced=0--> histogram has non-reduced time"
         write(0,*) "        =1--> histogram has reduced time"
         write(0,*) "           if data seems not so, stop is made"
         write(0,*) " actual number of time fraction data to be used"
         write(0,*) " nfraca=2-->up to T5, T10% is obtained "
         write(0,*)
     *         " nfraca=11-->up to T5,T10,...T90, T95% is obtained "
         write(0,*) " fitting is normally done by T10% so 2 is enough"

         write(0,*)
     *   "webonly=0-->input is made from mkLDD"
     *   "  (core region data exists) "
         write(0,*)
     *   "webonly=1-->input is made from FDD rawdata"
     *   "(no core region data) "
         write(0,*) " smooth      >0 max  # of smoothing"
         write(0,*) "          for LDD 100~500 may be ok"
         write(0,*) "          for FDD 3   may be ok"
         write(0,*) " ***** avoid smooth=0**** "
         write(0,*)
     *    " filename: path to the input -t.hist file(<256 chars)"
         stop 132                                        
      endif                                                
      call getarg(1, buf, status)
      read(buf,*) binorascii
      if(binorascii .ne. 1 .and. binorascii.ne. 2) then
         write(0,*) ' error input to binOrascii=',binorascii
         stop  444
      endif
      call getarg(2, buf, status)
      read(buf,*) reduced
      if(reduced .ne. 0 .and. reduced .ne. 1) then
         write(0,*) 'procTime: error input to reduced=',reduced
         stop  444
      endif

      call getarg(3, buf, status)
      read(buf,*) nfraca
      if(nfraca .lt. 1 .and. nfraca .gt. 11) then
         write(0,*) 'procTime: error input to nfraca=',nfraca
         stop  444
      endif

!                (defined in Zprivatesub) stdout 
      if(binorascii .eq. 1) then
!             ascii file
         call procTimeAscii(reduced, nfraca)
      else
!             bin file
         call procTimeBin(reduced, nfraca)
      endif
      end

      subroutine procTimeBin(reduced, nfraca)
      implicit none
#include "ZprivateSub.h"
#include "../../Hist/Z90histc.h"
#include "../../Hist/Z90histo.h"
#include "../../Hist/Z90hist1.h"
      character*80 buf
      integer reduced ! 0-->time is not reduced time
                      ! 1--> reduced time
      integer nfraca ! 2 --> T10% is used. (1~11)
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
      integer nbinhisto,  status, webonly
      integer kwhistReadAscii
      integer i
!      now  mode must be always 2 

      call getarg(4, buf, status)
      read(buf, *) webonly
      call getarg(5, buf, status)
      read(buf, *) smooth
      call getarg(6, filename, status)                          
!      write(0,*) status, filename

      open(fnoT, file=filename,
     *         iostat=status, access='sequential',
     *         form='unformatted', action='read')
      if(status .ne. 0 ) then
         write(0,*) ' cannot open file ='
         write(0,*) filename
         stop
      endif
      
!        skip time data for core region
      if(webonly .eq. 0) then
         do i = 1, ansites
            do code = 1, 4
               read( fnoT, end=100 ) histid0
               if( histid0 .ne. '#hist1' ) then
                  write(0,*) ' histogram is not 1D: ',histid0
                  stop
               endif
               call kwhistr(h10, fnoT, icon)
               call kwhists(h10, normf)
!                   deallocate
               call kwhistd( h10 )
            enddo
         enddo
      endif
!          web sector region
      do i = 1, ansites
         do code = 1, 4
            do ifai= 1, nfai
               do ir= 1, nrbin
                  read( fnoT, end=100 ) histid0
                  if( histid0 .ne. '#hist1' ) then
                     write(0,*) ' histogram is not 1-D: ',histid0
                     stop
                  endif
                  call kwhistr(h10, fnoT, icon)
!                     statistical calculation
                  call kwhists(h10, normf)
!                     time analysis: 6 is for stdout output
                  call procTime(h10, 6, nfraca, reduced, 
     *                          ir, ifai, code, i)
!                        deallocate
                  call kwhistd( h10 )
               enddo
            enddo  
         enddo   
      enddo
      close(fnoT)
      return
 100  continue
      write(0,*) ' unexpected EOF ' 
      stop 2345
!     **********************
      entry  procTimeAscii(reduced, nfraca)
!     **********************
      call getarg(4, buf, status)
      read(buf, *)  webonly
      call getarg(5, buf, status)
      read(buf, *) smooth
      call getarg(6, filename, status)                          
!      write(0,*) status, filename

      open(fnoT, file=filename,
     *     iostat=status, access='sequential',
     *     form='formatted', action='read')
      if(status .ne. 0 ) then
         write(0,*) ' cannot open file ='
         write(0,*) filename
         stop
      endif

      if(webonly .eq. 0 ) then
!        skip time data for core region
         do i = 1, ansites
            do code = 1, 4
               nbinhisto=kwhistReadAscii(h10, fnoT)
               if(nbinhisto .le. 0 ) then
                  write(0,*) ' ascii read  failed'
                  stop
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
                     stop
                  endif
!                     statistical calculation
                  call kwhists(h10, normf)
!                     time analysis
!                                6:  output on stdout
                  call procTime(h10, 6,  nfraca, reduced, 
     *                          ir, ifai, code, i)
!                        deallocate
                  call kwhistd( h10 )
               enddo
            enddo  
         enddo   
      enddo
      close(fnoT)
      end
