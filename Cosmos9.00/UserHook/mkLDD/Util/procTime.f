!
!
!     Suppose that the  normalzied integral (from left) time distribution
!     in a given web sector is like below  
!  1.0   ----------------------------------------*---- 
!                                          *
!                                      *
!  frac                         <-- *
!                                 * |
!                               *   |
!                             *     |
!                           *       |
!                         *         |
!                       *           |
!                     *             |
!                   *             
!                 *            
!                *
!              *                    tf           
! 0-------- *------------------------------------------> t
!    
!      We get set of 'tf' for given set of 'frac' which are
!      e.g, 0.05 0.1, 0.2,... 0.9, 0.95  (see fig). 
!
!      real*8 frac(nfrac)
!      real*8 frac/0.05d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0,
!     *                  0.7d0, 0.8d0, 0.9d0, 0.95d0/ 
!      real*8 tf(nfrac)  
!      ( nfrac=11 in this case. )
!
!      We get such tf's for every web sector and store
!      in an array tfary0; 
!      For a web sector with a given ir( lateral index) value 
!      we get first nfrac tf values. Then, we get similar one
!      for next web sector at next it.  So we get nrbin 'tf'
!      values for each 'frac'.  These are written to the
!      stdout with a header of "fai code layer" index.
!      The last data is indicated by "0 0  0" header index.
!
!      
!      
      subroutine procTime(h, fnotf, 
     *    nfraca, reduced, idxr, idxf, code, layer)
!         This treats one 1-D histogram with web index (idxr, idxf)
!         
      implicit none
#include "ZtimeAna.h"

#include "../../Hist/Z90histc.h"
#include "../../Hist/Z90histo.h"
#include "../../Hist/Z90hist1.h"
      
      type(histogram1):: h  ! input 1 D histogram
      integer fnotf  ! input file number for tf data  
      integer nfraca ! input.
                    !  1 to 11;  upto what % data is obtained.
                    !  For final fitting, T10% is ok so use 2.
                    !  1->T%5 2->T10%  9->T80% 10->T90% 11->T95%
                  ! 2 vs 11: 2 is 15% faster than 11. 
      integer reduced ! 0--> time is non-reduced time
                      ! 1--> reduced time
      integer idxr  ! web r bin index (1~nrbin=42)
      integer idxf  ! web fai bin index (1~nf) 
      integer code  ! ptcl code
      integer layer ! at which layer
      integer kwhistIxy
      integer maxsize  ! max histogram size
      integer n        ! actual histogram size
      
      integer i, j
      parameter (maxsize=3000)
      real*8 x(maxsize), y(maxsize)
      real*8 cgap
      integer icon
      integer idxr1, idxr2
      integer mode
      save
      mode=2  ! only this can be used

!       get (normalized) integral dist.
      n = kwhistIxy(h, x, y, maxsize)
      if(n .gt. maxsize) then
         write(0, *) ' too large histogram size=',n
         stop 1111
      endif
      if( idxr .eq. 1 ) then
         idxr1 = 0
         idxr2 = 0
      endif
      if(reduced .eq. 1 .and. idxf .eq. 1) then
!           check time<0 for fai=0 deg.
         if(n .gt. 0) then
            if(x(1) .lt. 0. ) then
               write(0,*) 
     *          'Although you gave reducedT="yes" in baseInfo'
               write(0,*)
     *           'I suspect that histgoram data has non reduced time'
               write(0,*) 'so I stop here'
               stop 999
            endif
         endif   
      endif

      if(smooth .gt. 0) then
!        smoothin parameter
         cgap = 0.1/(idxr/30.)**4
!          smoothing.( at least one trial is done.)
!          at large distances (idxr>30,  10 or more
!          trials) 
      
         call ksmooth(x, 1, y, 1,  n, 0, smooth, cgap, icon)
      else
         icon = 1
      endif

      if(icon .gt. 0 ) then
         if(idxr1 .eq. 0) idxr1=idxr
         idxr2=idxr
!            get 'tf' values corresponding to 'frac'
         call procTimeGetTf(x, y, n, frac, tf, nfraca)
         do i = 1, nfraca
            tfary0(idxr,  i) = tf(i)
         enddo
      endif
      if(idxr .eq. nrbin) then
!           all data at given fai has been obtained.
!           data is between idxr1 to idxr2 for r
!           do fitting        and save the result in tfary 

         write(fnotf,'(3i4)')  idxf, code, layer
         if(mode .eq. 1) then
            do i = idxr1, idxr2
               write(fnotf,'(i3,1p11E11.3)') 
     *         i, (tfary(i, j), j=1, nfraca)
            enddo
            write(fnotf,
     *       '("0 0 0 0 0 0 0 0 0 0 0 0 ")')
         elseif(mode .eq. 2) then
            do i = idxr1, idxr2
               write(fnotf,'(i3,  1p11E11.3)')
     *            i, (tfary0(i, j), j=1, nfraca)
            enddo
            write(fnotf,
     *       '("0 0 0 0 0 0 0 0 0 0 0 0 ")')
         elseif(mode .eq. 3) then
            do i = idxr1, idxr2
               write(fnotf,'(i3,1p22E11.3)')
     *          i, (tfary(i, j), j=1, nfraca),
     *           (tfary0(i,j), j=1, nfraca)
            enddo
            write(fnotf,
     *       '("0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0")')
         endif
      endif
      end

      subroutine procTimeGetTf(x, y, n, frac, tf, nfrac)
      implicit none
      integer n
      real*8 x(n)
      real*8 y(n)
      integer nfrac
      real*8 frac(nfrac)
      real*8 tf(nfrac)

      real*8 error
      integer i, j
      integer np
      parameter (np = 3)  ! use np points for interpolation

      do i = 1,  nfrac
         call kpolintpFE(y, 1, x, 1, n,  np,
     *      frac(i),  tf(i), error)     
      enddo
      end
