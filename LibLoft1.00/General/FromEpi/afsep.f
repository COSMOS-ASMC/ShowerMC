       subroutine afsep(io)
       implicit none
       integer io
          character* 10  sep
!             *** until loop*** 
             do while (.true.)
               read(io, '(a)') sep
             if         (sep .eq. '----------')
     *                          goto 10
             enddo
   10        continue
       end
