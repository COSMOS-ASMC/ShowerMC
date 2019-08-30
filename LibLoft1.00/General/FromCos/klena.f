!      ***************************************************************
!      *                                                             *
!      * klena: actual length of character string
!      *                                                             *
!      ********************** tested. 87.01.16*******************k.k**
!
!       /usage/   l=klena(cha)
!         e.g.    character*10 cha
!                 cha='ab c d'
!                 call sub(cha)
!                 ....
!                 subroutine sub(cha)
!                 character*(*) cha
!                 l=klena(cha)
!                   l will be 6
!    note:  if no character dec. is given for string, klena=0 will
!           result.  e.g.  a='abc',  klena(a). but klena('abc') is ok.
!
       integer function klena(cha)
       implicit none
       character*(*) cha
       integer l

       l=len(cha)

       if(l <=  0 .or. cha == " " ) then
          klena=0
       else
          do   while (cha(l:l) .eq. ' ')
             l=l-1
          enddo
          klena=l
       endif
       end
