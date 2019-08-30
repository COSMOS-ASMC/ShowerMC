!       ***************************************************************
!       *
!       *  epExpot:  compute excitation potential (in eV)
!       *            consisting of n different atoms.
!          
!         ---------------------------------
!          must be called after epGetEffZA
!         ---------------------------------
!       *
!
! /usage/ call epExpot(media)
!      
!           see p.r. b vol3 (1971)3681 sternheimer.  c

       subroutine epExpot(media)
       implicit none
#include "Zmedia.h"

       type(epmedia):: media  !  input.  data  given by basic table must
                            !          be ready
                            !  output. media.I. ionization potential in GeV

       integer n

       integer i, j
       real*8 sumli, expi, temp

       real expot(92)
!
!   From  http://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html
!    expot(Z) is excitation poential for Z in eV. Sep. 7, 2002
!   These are the same as Seltzer & Berger(see below) and also
!   Groom et al. Atomic data and Nucl. data Tables vol.78, p.183 (2001)
!      Some are for gas. (H,He etc).  The error cited in Seltzer & Berger
!      is typically a few percent.
       data expot/
     1 19.2,  41.8,  40.0,  63.7,  76.0,   78.0,  82.0,  95.0, 115.0,
     2 137.0, 149.0, 156.0, 166.0, 173.0, 173.0, 180.0, 174.0, 188.0,
     3 190.0, 191.0, 216.0, 233.0, 245.0, 257.0, 272.0, 286.0, 297.0,
     4 311.0, 322.0, 330.0, 334.0, 350.0, 347.0, 348.0, 343.0, 352.0,
     5 363.0, 366.0, 379.0, 393.0, 417.0, 424.0, 428.0, 441.0, 449.0,
     6 470.0, 470.0, 469.0, 488.0, 488.0, 487.0, 485.0, 491.0, 482.0, 
     7 488.0, 491.0, 501.0, 523.0, 535.0, 546.0, 560.0, 574.0, 580.0,
     8 591.0, 614.0, 628.0, 650.0, 658.0, 674.0, 684.0, 694.0, 705.0,
     9 718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0, 800.0, 810.0,
     a 823.0, 823.0, 830.0, 825.0, 794.0, 827.0, 826.0, 841.0, 847.0,
     b 878.0, 890.0/


       sumli=0.



       n = media%noOfElem

       do  i = 1, n
!          if(media.elem(i).Z .eq. 1) then
!             expi=18.7
!          elseif(media.elem(i).Z .lt. 13.0) then
!             expi=13.*media.elem(i).Z
!          else
!             expi=media.elem(i).Z *( 9.76 + 58.8*
!     *            media.elem(i).Z**(-1.19))
!          endif
!
          j  =  media%elem(i)%Z 
          expi = expot( j )   ! OSX absoft must use j for inteface
          if(n .gt. 1) then
!            we use 13 % rule ; though good only for condensed media
!            ref: Seltzer & Berger, Int.J. Rad. Isot. vol.33.p.1180, 1982 
             expi = expi*1.13  ! from v8.80
          endif
!                 bef. 8.80; but without /media.ZbyAeff after enddo
!          sumli=sumli +      
!     *       media.elem(i).Z* media.No(i)/media.Z * log(expi)
!                 Zi * Ni/(sum Ni*Zi) *log(expi)
!                from 8.80
          sumli=sumli + 
     *           media%w(i)*media%elem(i)%Z/media%elem(i)%A*log(expi)
       enddo
       sumli = sumli /media%ZbyAeff
       temp = exp(sumli) * 1.d-9  ! to GeV

       if( media%format .ne. 2 .or. media%I .eq. -100.) then
          media%I = temp
       else
!            'I' sholud  have been given by input
          if( abs(media%I-temp)/media%I .gt. 0.01) then
!               issue message only if diff > 20%
             write(0,*) ' For media=', media%name
             write(0,*) ' excitation potential given by input =',media%I
             write(0,*) ' computed value inside is =', temp, ' eV'
             if( abs(media%I-temp)/media%I .gt. 0.5)  then
                write(0,*) ' diff is too large; check the input'
                stop
             endif
             write(0,*) ' we use input value'
          endif
       endif
       end
