c
c     test cthik2len.f
c
c  
      implicit none
      real*8  len, z, cosz, thick, cutt
      real*8  cvh2thick, tempz, clen2thickEx
      integer ios, jcut

      do while(.true.)
         z = 0.d3
         cosz = 1.d0
         thick = 1000.d0
         write(0,*) " enter  z, cosz, thick", z, cosz, thick
         read(*, *, iostat=ios) z, cosz, thick
         if(ios .ne. 0) stop 
         call cthick2len(z, cosz, thick, len, cutt, jcut)
         write(*, *) " cutt= ", cutt, " jcut=", jcut
c
         if(abs(cosz)  .eq. 1.) then
            write(*,*) (cvh2thick(z-len*cosz) - cvh2thick(z) )* cosz
            tempz = z - len*cosz
            write(*, *) " new length = ", len, " height =", tempz
            write(*, *) " thickness at z-len*cos=", cvh2thick(tempz)
            write(*,*)" thickness=",
     *      (  cvh2thick(tempz) - cvh2thick(z) )* cosz
         else
            thick= clen2thickEx(z, cosz, len)
            write(*, *) ' thickness for len=', thick
         endif
      enddo
      end
