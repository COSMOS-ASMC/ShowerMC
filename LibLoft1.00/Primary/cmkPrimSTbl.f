! *******************************************************
      subroutine cmkPrimSTbl(each, icon)
!         make primary energy sampling table for a given
!         component.
!     *******************************************************
      implicit none

#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryc.h"
      type(component)::  each
      integer  icon
!
      integer k, nseg, l

      character*150 msg
      integer firstseg, i, lastseg
      real*8 cut, cut2
!
      nseg = each%no_of_seg    ! # of segmets = datapoints -1
      cut = each%cut
      cut2 =each%cut2
!         first convert flux into normal expression
!         energy unit is kept as input
!         original each flux = (dI/dE)*E**flatterer
!
!        check energy order.         
      do l = 1, nseg
         if(each%energy(l) .ge. each%energy(l+1)) then
             write(*,*) each%label,'-th component=',each%symb,
     *       ' energy is not in ascending order.'
             write(*,*) ' It is in the ', l, '-th segment'
          endif
      enddo
!         flux without multiplication of E**flatterer
      do l = 1, nseg+1
         each%flux(l) = each%flux(l)/each%energy(l)**each%flatterer
      enddo
      firstseg = 1
      lastseg = nseg

!           now flux = dI/dE = const *E**(-beta). get beta next
      if(each%diff_or_inte .eq. 'd') then
!             get beta in E**(-beta)dE in each segment
           do k =1, nseg
              each%beta(k) = - log(each%flux(k)/each%flux(k+1))
     *                        /
     *                        log(each%energy(k)/each%energy(k+1))
           enddo
!               assume last segment is the end of spectrum
           if(nseg .ge. 1) each%beta(nseg+1) = 100.
!
!          -----------------
!          find segment i <= cut < i+1, if cut > 0
!          and adjust each.emin.
!
        if(cut .gt. 0. .and.  nseg .ge. 1 ) then
           if(cut .ge. each%energy(nseg+1)) then
              call cerrorMsg(
     *        'the primary minmum cut >= table max',0)
           endif
           do i = nseg, 1, -1
              if(cut .ge. each%energy(i)) then
                 firstseg = i
                 each%emin = cut
                 goto 100
              endif
           enddo
           each%emin = each%energy(1)
           firstseg = 1
!            nothing to do
           goto 200
 100       continue
!             flux  at cut
           k = firstseg
           each%flux(k) = each%flux(k)* 
     *    ( cut/each%energy(k))**(-each%beta(k))
           each%energy(k) = cut
        endif
 200    continue
!          ------------------
!          find segment i < cut2 <= i+1, if cut2 > 0
!          and adjust each.emax,
!


        if(cut2 .gt. 0. .and.  nseg .ge. 1 ) then
           if(cut2 .le. each%energy(1)) then
              call cerrorMsg(
     *        'the primary max cut <= table min',0)
           endif
           do i =  1, nseg
              if(cut2 .le. each%energy(i+1)) then
                 lastseg = i
                 each%emax=  cut2
                 goto 300
              endif
           enddo
           each%emax = each%energy(nseg+1)
           lastseg = nseg
!            nothing to do
           goto 400
 300       continue
!             flux  at cut2
           k = lastseg
           each%flux(k+1) = each%flux(k)* 
     *    ( cut2/each%energy(k))**(-each%beta(k))
           each%energy(k+1) = cut2
        endif
 400    continue

!              get integral distribution
        if(nseg .ge. 1) then
!           assume spectrum is cut at each.emax
           each%norm_inte(lastseg+1)= 0.
           do k = nseg+1, lastseg+1, -1
              each%norm_inte(k) = 0.
           enddo
           do k = lastseg, firstseg, -1
              if(abs(each%beta(k)-1.d0) .gt. 1.d-6) then
                   each%norm_inte(k) = each%norm_inte(k+1) +
     *             (  each%energy(k+1)*each%flux(k+1)
     *            -each%energy(k)*each%flux(k) )/(1.d0-each%beta(k))
              else
                  each%norm_inte(k) =  each%norm_inte(k+1) +
     *                 each%flux(k)* each%energy(k)*
     *                 log(each%energy(k+1)/each%energy(k))
              endif
           enddo
        endif
      elseif(each%diff_or_inte .eq. 'i') then
!               integral. first make power of integral spectrum
         do k =1, nseg
            each%beta(k) = - log(each%flux(k)/
     *           each%flux(k+1))
     *           /
     *           log(each%energy(k)/each%energy(k+1))
         enddo
         if(nseg .ge. 1) then
            each%beta(nseg + 1) = each%beta(nseg)
         endif
!               move flux into norm_inte
         do k = 1, nseg+1
            each%norm_inte(k) = each%flux(k)
         enddo
      else
         write(*,*) ' error specification of diff/integral=',
     *        each%diff_or_inte
         icon = 1
         goto 900
      endif

!                total integral value
      if(nseg .ge. 1) then
         each%inte_value = each%norm_inte(firstseg)
      else
!                  delta funcition
         each%inte_value = each%flux(1)
      endif
!               normalize the flux
      do k = firstseg, lastseg
         each%norm_inte(k) = each%norm_inte(k)/each%inte_value
      enddo
 900  continue
      end





