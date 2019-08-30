      subroutine kxplbxNew(x0, y0, z0,  l, m, n,  a, b,  c,
     *            el,  icon)
!          if kxplbx defined in KKlib has some problems
!          use this one (rename is of course needed)
!       kxplbx defined in KKlib is little bit faster
!       (~7%).
!
      implicit none
      real(8),intent(in):: x0, y0, z0
      real(8),intent(in):: l, m, n
      real(8),intent(in):: a, b, c
      real(8),intent(out):: el
      integer,intent(out):: icon
!
      real(8)::p1(3), p2(3), p3(3), p4(3)
!
      real(8):: sp(3), dir(3)
      logical::insideflag
      integer::icon1, icon2
      real(8):: leng, lenga(2)
      logical kInsideBox ! external func.
      integer,parameter::face(6)=(/1,6,2,5,3,4/)
      integer::i,  nx
      real(8)::dummy

      real(8),parameter:: eps=-1.d-10

      insideFlag = kInsideBox(x0, y0, z0, a, b, c)
      nx = 0
      sp(1) = x0
      sp(2) = y0
      sp(3) = z0
      dir(1) = l
      dir(2) = m
      dir(3) = n
      
      do i = 1, 6
         call kFormBoxFace( face(i), a,b,c, p1, p2, p3, p4) 
         call epxpLand4vp(p1, p2, p3, p4, sp, dir,
     *                     leng, icon1, icon2)
         if( icon1 <= 4 ) then
            if( insideFlag ) then
               if( leng >  0.d0 ) then
                  nx = nx +1
                  lenga(nx)= leng
                  exit
               endif
            else
               if( leng > 0.d0 ) then
                  nx = nx + 1 
                  lenga(nx) = leng
                  if( nx == 2  ) exit
               endif
            endif
         endif
      enddo
      if( nx == 0 ) then
         if( insideFlag ) then
            write(0,*) 'in kxplbx:  strange inside but no xp' 
            write(0,*) ' x0, y0, z0=',x0, y0, z0
            write(0,*) ' l, m, n=', l, m, n
            write(0,*) ' a, b,  c=', a, b,  c
            stop
         endif
         icon = -1
      elseif( nx == 1 ) then
         if( insideFlag )  then
            icon = 0
            el = lenga(1)
         else
            ! for safety message
            write(0,*) ' strange in kxplbx'
            call epfordebug('kxplbx')
            call cbackTrace(dummy)
            stop
         endif
      else
         if( insideFlag ) then
            write(0,*)' strange inside but xp=2'
            write(0,*) ' leng =', lenga(1:nx)
            call epfordebug('kxplbx')
            call cbackTrace(dummy)
            stop
            el = maxval(lenga(:))
            icon = 0
         else
            el = minval(lenga(1:nx))
            if(el < 0.d0) then
               write(0,*) ' strange in kxplbx; 0-->el '
               write(0,*) ' lenga=',lenga,' nx= ', nx
               write(0,*) ' a,b,c =', a, b, c
               write(0,*) ' sp =', sp
               write(0,*) ' dir =', dir
               call epfordebug('kxplbx')
               call cbackTrace(dummy)
               stop
            endif
            icon = 1
         endif
      endif
      end  subroutine kxplbxNew
