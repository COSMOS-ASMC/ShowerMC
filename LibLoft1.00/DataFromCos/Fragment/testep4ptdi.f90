      implicit none
      real(8)::xa(2)=(/0., 1./)
      real(8)::ya(2)=(/3., 0.5/)
!      real(8)::f(2,2)=(/1., 2.5, 2., 3./)
      real(8)::f(2,2)=(/1., 0., 2., 3./)
      real(8):: hx=1.d0
      real(8):: hy=-2.5d0
      real(8):: x0= 0.
      real(8):: y0= 3.
      real(8):: ans, x, y 
      
!      x = 0.5d0
!      y = 2.0d0
!      call ep4ptdi(f,   x0, y0, hx, hy, x, y, ans)
!      write(*,*) x, y, ans
!      x = 0.9d0
!      y = 2.8d0
!      call ep4ptdi(f,   x0, y0, hx, hy, x, y, ans)
!      write(*,*) x, y, ans
!      x = 0.2
!      y = 2.5
!      call ep4ptdi(f,   x0, y0, hx, hy, x, y, ans)
!      write(*,*) x, y, ans
!      x = 0.95
!      y = 0.6
!      call ep4ptdi(f,   x0, y0, hx, hy, x, y, ans)
!      write(*,*) x, y, ans
      x = 0.2
      y = 1.0
      call ep4ptdi(f,   x0, y0, hx, hy, x, y, ans)
      write(*,*) x, y, ans
      end
!     ****************************************************************
!     *                                                              *
!     * ep4ptdi:  4-point two dimensional interpolation              *
!     *   special for epSampNHR
!     *          each of two points maybe the same points
!     *   if function vlaue is 0, not used.
!     ****************************************************************
!
!   /usage/
!
!       call ep4ptdi(f, x0, y0, hx, hy, x, y, ans)
!
!     f containes the function values at (x0,y0), (x0+hx, y0+hy)
!     ans gets the value of the funtion at (x,y)
!     linear interpolation is tried.
!     hx can be 0, if x is the same as x0
!     hy can be 0, if y is the same as y0
!     if f is 0, only non zero value is used.
!     if all f is 0, ans=0.
      subroutine ep4ptdi(f,  x0, y0, hx, hy, x, y, ans)
      implicit none
      real*8 f(2,2),  x0, y0, hx, hy, x, y, ans
!
      integer i, j, i1, j1
      real(8):: a, b, pL, pH, q, pL1, pH1, q1
      real(8):: fij, fi1j, fij1, fi1j1, fa, fb

      if(hx /= 0.0) then
         a=(x-x0)/hx
      else
         a = 0.
         if(x /= x0) then
            write(0,*) ' x=',x, ' != x0=',x0,' while hx=0'
            write(0,*) ' in ep4ptdi'
            stop
         endif
      endif
      if( hy /= 0.0) then
         b = (y-y0)/hy
      else
         b = 0.
         if(y /= y0) then
            write(0,*) ' y=',y, ' != y0=',y0,' while hy=0'
            write(0,*) ' in ep4ptdi'
            stop
         endif
      endif

      i=1
      j=1
      pL=a
      pH=a  
      pL1=1.-pL
      pH1=1.-pH
      q=b
      q1=1.-q
      if( hx == 0.) then
         i1 = i
      else
         i1 = i + 1
      endif
      if( hy  == 0.) then
         j1 = j
      else
         j1 = j + 1
      endif
      fij = f(i,j)
      fi1j= f(i1,j)
      fij1= f(i,j1)
      fi1j1= f(i1,j1)  

      write(0,*) ' fij:',fij, 'fi1j=',fi1j
      write(0,*) ' fij1:',fij1, 'fi1j1=',fi1j1
      if( fij == 0.) then
         pL = 1.
         pL1 =0.
      endif
      if( fi1j == 0.) then
         pL1 = 1
         pL =  0.
      endif
      if( fij1 == 0. ) then
         pH = 1.
         pH1 = 0.
      endif
      if( fi1j1 == 0.) then
         pH = 0.
         pH1 = 1.
      endif
      write(0,*) ' pL, pL1 =', pL, pL1
      write(0,*) ' pH, pH1 =', pH, pH1
      fa =  fij*pL1 + fi1j*pL
      fb =  fij1*pH1 + fi1j1*pH 
      write(*,*) x,y0, fa
      write(*,*) x,y0+hy, fb

      if( fa == 0.) then
         q1 =0.
         q = 1.
      elseif(fb == 0.) then
         q1 = 1.
         q = 0.
      endif
      write(0, *) 'q=',q, 'q1=',q1
      write(0,*) ' fa=',fa, 'fb =',fb

      ans=  fa  * q1 +  fb* q
   end subroutine ep4ptdi

