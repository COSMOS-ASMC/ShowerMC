!     
!            weighted histograming without dynamic memory allocation
!      Usage:  kwhisti:   initialization for one histogram
!              kwhistc:   clear histogram area
!              kwhist:    take histogram
!              kwhists:   compute statistical result 
!              kwhistp:   print statistical resutl
!
!
      subroutine kwhisti(h, ixmin, ibin,  dummy, itklg )
      implicit none
!         initialize 
      integer dummy   !  not used. for the compat. with f90 version.
      real ixmin     ! input. xmin. not in log even itklg ==1
      real ibin      ! input. bin.  If log, bin is for log10
      integer itklg    ! input.  bit 1: 0--> not take log10
                       !                1--> take log10
                       !             2: 0--> ixmin is the min of bin
                       !                1--> ixmin is the center of the min bin
                       !             3: 0--> neglect underflow   
                       !                1--> underflow is put in min bin
                       !                     mean bin value is affected by
                       !                     those with underflowed values
                       !             4: 0--> neglect overflow
                       !                1--> overflow is put in the highest bin
                       !                     mean bin value is affected by
                       !                     overflowed ones
!     ******************
#include "Zhist.f"
      type(histogram):: h
!     ====================      
      real inorm   !  input. used in the normalization as dN/dx/inorm
                   !  if this is 0,  area normalization is tried.
!
      real  x, w 
      real  xx
      integer isum,  i, ndiv
      logical asmax
      real isumw
      real dx
      character*(*)  id
      
      h.tklg = ( itklg-(itklg/2)*2  ) .ne. 0
      h.cent = ( (itklg/2)*2 - (itklg/4)*4)/2
      h.ufl  = ( (itklg/4)*4 - (itklg/8)*8 ) .ne. 0
      h.ofl  = (  (itklg/8)*8 - (itklg/16)*16 ) .ne. 0

      h.xmin = ixmin  ! (not used at present)
      asmax =   (  (itklg/16)*16 - (itklg/32)*32 ) .ne. 0
      if(asmax) then
         if(ixmin .ge. ibin) then
            write(0,*) ' ibin is regarded as ixmax but <= ixmin'
            stop 99999
         else
            if( h.cent .eq. 1 ) then
               ndiv= ibin - 1
            else
               ndiv = ibin
            endif
            if(h.tklg) then
               h.bin = log10(ibin/ixmin)/ndiv
            else
               h.bin = (ibin - ixmin )/ndiv
            endif
         endif
      else
         h.bin = ibin
      endif
      

      h.nhist = nhistp
      if( h.tklg ) then
         if(ixmin <= 0.0)  then
            write(0,
     *       '("min must be > 0 for log option")')
            stop
         endif
         h.xm = log10(ixmin) - h.cent * h.bin/2.
         h.inc = 10.**h.bin
      else 
         h.xm = ixmin - h.cent *  h.bin/2.
         h.inc  = h.bin
      endif
      return

!     *****************
      entry kwhistc(h)
!     *****************
      do i = 1, h.nhist
         h.dn(i) = 0.
         h.dnw(i) = 0.
         h.mean(i) = 0.
      enddo
      return

!    *************************
      entry kwhist( h, x, w )
!    *************************
      if( h.tklg  .and. x .le. 0.) then
!         neglect this data
      else
         if( h.tklg  ) then
            xx = log10(x)
         else
            xx = x
         endif
         i = ( xx-h.xm ) / h.bin  + 1
         if(i .le. 0 .and. h.ufl ) then
            i = 1
         elseif(i .gt. h.nhist .and. h.ofl ) then
            i = h.nhist
         endif
         if(i .ge. 1 .and.  i  .le. h.nhist )  then
            h.dn(i) = h.dn(i)  + 1
            h.dnw(i) = h.dnw(i) + w
            h.mean(i) = h.mean(i) + x*w
         endif
      endif
      return

!     ***********************
      entry kwhists( h, inorm )
!     ************* take statistics
      h.norm = inorm
      h.imin = 1
      
      do while(h.imin .lt. h.nhist .and. h.dn(h.imin) .eq. 0) 
         h.imin = h.imin + 1
      enddo

      h.imax = h.nhist
      do while ( h.imax .gt. 1 .and.  h.dn(h.imax) .eq.  0)  
         h.imax = h.imax -1
      enddo

      h.sum= 0
      do i = h.imin, h.imax
         h.sum =  h.sum + h.dn(i)
      enddo

      h.sumw = 0.
      do i = h.imin, h.imax
         h.sumw = h.sumw +  h.dnw(i)
      enddo
!           bin center value
      if(h.tklg ) then
         xx = 10.**( h.xm + h.bin/2.) * h.inc**(h.imin-1)
      else
         xx = h.xm + h.bin/2. + h.inc*(h.imin-1)
      endif

      dx = h.bin      

      if(h.norm .eq. 0. .and. h.sumw .gt. 0.) then
         h.norm =  h.sumw
      else(h.norm .eq. 0. ) then
         h.norm = 1.
      endif

      do i = h.imin, h.imax
         if(h.tklg ) then
            dx  = 10.0**(h.xm + i * h.bin) - 10.0**(h.xm + (i-1)*h.bin)
         endif
         if(h.dnw(i) .eq. 0) then
            h.mean(i) = xx
         else
            h.mean(i) = h.mean(i)/h.dnw(i)
         endif
         h.dndx(i) = h.dnw(i)/dx/h.norm
         if(h.tklg ) then
            xx = xx * h.inc
         else
            xx = xx + h.inc
         endif
      enddo
      return
!     *********************
      entry kwhistp(h, id)
!     ****************print  hist
      isum = h.sum
      isumw = h.sumw


      if(h.tklg ) then
         xx =10.0**( h.xm + h.bin/2) * h.inc**(h.imin-1)
      else
         xx = (h.xm +  h.bin/2) + h.inc*(h.imin-1)
      endif

      do i = h.imin, h.imax
         write(*, '(a, " ", i5,  4g14.4, i10, i10, g14.4)') id, i,
     *   xx, h.dndx(i), h.dnw(i), isumw, h.dn(i), isum,  h.mean(i)
         isum = isum - h.dn(i)
         isumw = isumw -  h.dnw(i)
         if( h.tklg ) then
            xx =  xx * h.inc 
         else
            xx = xx + h.inc
         endif
      enddo
      end






