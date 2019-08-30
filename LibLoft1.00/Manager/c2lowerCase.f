!c         test c2lowerCase
!      character*15  cu/'Ab~} C-_De'/
!      character*15  cl/'Xc~+D-_DE'/
!       write(*, *) cu, cl, '/'
!       call c2lowerCase(cu, cl)
!       write(*, *) cu, cl, '/'
!       call c2upperCase(cl, cu)
!       write(*, *) cl, cu, '/'
!       end
      subroutine c2lowerCase(cu, cl)
!       
!        convert [A-Z] in "cu" into [a-z] and put it in "cl"
!        This is for standard ASCII.
!
!    cu: character string. input
!    cl: //              . output.  cl may be cu.
!
!             if cl is shorter than cu, overflowing part
!          is not moved.
!
      implicit none
      character*(*) cu, cl
!
      integer i, ic, imax
!

      imax =  min(len(cu), len(cl))
      do i = 1, imax
         ic = ichar( cu(i:i) )
         if(  ic .ge. 65 .and. ic .le. 90 ) then  ! A-Z
             cl(i:i) = char(ic+32)
         else
             cl(i:i) = cu(i:i)
         endif
      enddo
      end
      subroutine c2upperCase(cl, cu)
!       
!        convert [a-z] in "cl" into [A-Z] and put it in "cu"
!        This is for standard ASCII.
!
!    cl: character string. input
!    cu: //              . output.  cu may be cl.
!
!             if cu is shorter than cl, overflowing part
!          is not moved.
!
      implicit none
      character*(*) cl, cu
!
      integer i, ic, imax
!

      imax =  min(len(cu), len(cl))
      do i = 1, imax
         ic = ichar( cl(i:i) )
         if(  ic .ge. 97 .and. ic .le. 122 ) then  ! a-z
             cu(i:i) = char(ic-32)
         else
             cu(i:i) = cl(i:i)
         endif
      enddo
      end

