!      ************************************************************
!      *
!      *  kseblk: supress extra blanks
!      *          and replace specified character by blank
!      *
!      *********************** tested 87.01.29 ***************k.k**
!
!      /usage/  call kseblk(text, c, lc)
!
!    input   text:  character string
!               c:  one character to be replaced to blank
!    output
!            text:  character string.  first blank is supressed and
!                   then c is replaced to blank.
!              lc:  resultant text length
!
       subroutine kseblk(text, c, lc)
       implicit none
       integer lc
!
       character*(*) text
       character*1 c

       integer klena, l, i, lg
!
       l=klena(text)
       lc=0
        do   i=1, l
           if(text(i:i) .ne. ' ') then
               lc=lc+1
               if(text(i:i) .eq. c) then
                    text(lc:lc)=' '
               else
                    text(lc:lc)=text(i:i)
               endif
           endif
        enddo
       lg=len(text)
       if(lg .gt. lc) then
           text(lc+1:lg)=' '
       endif
       end
