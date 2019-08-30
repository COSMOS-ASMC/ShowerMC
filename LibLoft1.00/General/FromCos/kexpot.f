!       ***************************************************************
!       *
!       *   kexpot:  compute excitation potential ep (in eV)
!       *     consisting of n different atoms.
!       *
!       *   see P.R  B vol3 (1971)3681 sternheimer.  
!       *
!       ***************************** tested 87.09.19 *****************
!
! /usage/ call kexpot(za, aa, rela, n, zi, ai, ep)
!
!      Example for C3H8
!      parameter (n=2)
!      real*8 aa(n), za(n), rela(n), zi, zi, ep
!      data aa/1., 12./, za/1., 6./, rela/8., 3./
!      call kexpot(za, aa, rela, n, zi, ai, ep)
!      write(*,*) ' ep=', ep
!      end
!
        subroutine kexpot(za, aa, rela, n, zi, ai, ep)
         implicit none
         integer n
        real*8 za(n) ! input. za(i) is the charge of the i-th atom
        real*8 aa(n) ! input. aa(i) the mass # of i-th atom
        real*8 rela(n) ! input. rela(i) the relative number of i-th atom
         real*8 zi      ! output. weighted sum of za
        real*8 ai      ! output. weighted sum of aa
        real*8 ep      ! output. Excitention potential in eV

!
         integer i
         real*8 sumli, expi

        ai = 0.
        zi = 0.
        do   i = 1, n
             ai = ai + aa(i)*rela(i)
             zi = zi + za(i)*rela(i)
        enddo
!
        sumli=0.
        do   i=1, n
             if(za(i) .eq. 1) then
                   expi=18.7
             elseif(za(i) .lt. 13.) then
                   expi=13.*za(i)
             else
                   expi=za(i) *( 9.76 + 58.8*za(i)**(-1.19))
             endif
             sumli=sumli + za(i)*rela(i)/zi * log(expi)
        enddo
        ep=exp(sumli)
       end
