      implicit none
      include "Zprivate2.f"

      character*120  tracefile
      integer leng, icon, i, outtype
      real*4  mint, maxt, dt, t1, t2
      logical dothin, split
      integer klena

      mint = 000.
      maxt = 20000.
      dt = 10.
      pixel = dt
      read(*,*)  
     *  tracefile, dir,  mint, maxt, dt, pixel, dothin, split,
     *  outtype, mulalpha

      write(0,*)
     *  tracefile(1:klena(tracefile)), dir(1:klena(dir)),
     *  mint, maxt, dt, pixel, dothin, split,
     *  outtype, mulalpha

      call copenfw2(11, tracefile, 1, icon)
      if(icon .ne. 1) then
         write(0,*) tracefile, ' cannot be opened'
         stop 1111
      endif

      t2 = mint
      do while (t2 .lt. maxt)
         t1 = t2
         t2 = min(t1+ dt*n , maxt)
         jmax = 0
         jmin = 10000000
         call slice(t1, t2, dt)
         write(0,*) ' t1=',t1, 't2=',t2, ' slice ended'
         if(dothin) then
            call thinning(jmin, jmax)
            write(0,*) ' thinning ended'
         endif

         if(split) then
            if( ibits(outtype, 0, 1) .ne. 0) then
               call mkfiles
            endif
            if( ibits(outtype, 1, 1) .ne. 0) then
               call mkskelfiles
            endif
            write(0,*) ' all output finished'
         else
            if(ibits(outtype, 0, 1) .ne. 0) then
               call mkafile
            endif
            if(ibits(outtype, 1, 1) .ne. 0) then
               call mkaskelfile
            endif
            write(0,*) ' ouput to a file ended'
         endif
         rewind(11)
      enddo
      end
!     *********************************
!       read one trace file and  slice data by constant time surface.
!       data is accumulated in x(i,j)  etc.  where i is the index of 
!       number of ptcls in the given time plane. j the index of the
!       time plane.  
!     *********************************
      subroutine slice(mint, maxt, dt)
      implicit none
      include "Zprivate2.f"

      real*4  mint, maxt, dt
      integer i, j
      real*4  t,  k
      integer*2 code, chg 
      integer it1, it2
      real*4  xx1, yy1, zz1, time1
      real*4  xx2, yy2, zz2, time2

      character* 120 input

      do j = 1, n
         idx(j) = 0
      enddo
      
      do while (.true.)
         input = ' '
         read(11, '(a)', end=1000 )  input
         if( input .ne. " " ) then
            read(input, * ) xx1, yy1, zz1,  time1, code, chg
            do while( input .ne. " ") 
               read(11, '(a)', end =1000)  input
               if( input .ne. " " ) then
                  read(input, * ) xx2, yy2, zz2,  time2, code, chg
                  it1 = (time1- mint)/dt +1
                  it2 = (time2- mint)/dt 
                  if(it1 .lt. 1) goto 900
                  if(it2 .gt. n) goto 900
                  if(time2 .gt. maxt) goto 900
                  jmin = min(jmin, it1)
                  jmax = max(jmax, it2)
                  do  j = it1, it2
                     t = j*dt + mint
                     if( idx(j) .ge. maxp ) then
                        write(0,*)
     *                    'ptcls at time=', t , ' > ',  maxp
                        write(0,*)
     *                    'try to thin ptcls in the same pixel'
                        call thinning(j,j)
                     endif
                     idx(j) = idx(j) + 1
                     if(idx(j) .gt. maxp) then
                        write(0,*) 'ptcls > ',  maxp
                        write(0,*) 'thinning not worked'
                        stop 1234
                     endif
                        
!                           (x2-x1)/(t2-t1) *(t-t1) +x1
                     k = (t - time1) /(time2-time1)
                     x(idx(j), j) = (xx2-xx1) * k + xx1
                     y(idx(j), j) = (yy2-yy1) * k + yy1
                     z(idx(j), j) = (zz2-zz1) * k + zz1
                     codex(idx(j), j) = code
                     chgx(idx(j), j) = chg
                  enddo   
 900              continue
                  xx1 = xx2
                  yy1 = yy2
                  zz1 = zz2
                  time1 = time2
               endif
            enddo
         endif
      enddo
 1000 continue
      end
!     ***********************************
!      do thinning if two or more data fall in the dame pixel.
!     **********************************
      subroutine thinning(i1, i2)
      implicit none
      include "Zprivate2.f"
!        culling overlapping points in a pixel
!           this is very rough thinning
!     sorr each x
      integer i1, i2
      integer nc, j, i, k1, k2, l, nc0
      logical  dothin

      do j = i1, i2
         dothin = .false.
         if( idx(j) .gt. 0) then
            call kqsortr(x(1,j), si, idx(j))
!            call ksortinv(si, idx(j))
            do i = 1, idx(j)-1
               k1 = si(i)
               if(k1  .lt. 0) goto 20
               do l = i+1,  idx(j)
                  k2 = si(l)
                  if(k2 .lt. 0) goto 10
                  if(abs(x(k1, j)- x(k2,j)) .gt. pixel ) goto 20
!                         since xx is sorted,  we may safe to skip
!                         further check
!
                  if(abs(y(k1, j)- y(k2,j)) .gt. pixel ) goto 10
                  if(abs(z(k1, j)- z(k2,j)) .gt. pixel ) goto 10
!                   remove  k2
                  si(l) = - k2
                  dothin = .true.
 10               continue
               enddo
 20            continue
            enddo
            if(dothin) then
               nc0 = idx(j)
               call mvdata(x, j, nc)
               call mvdata(y, j, nc)
               call mvdata(z, j, nc)
               call mvdatai(codex, j, nc)
               call mvdatai(chgx, j, nc)  
               idx(j) = nc
!///////////// 
               write(0,*) ' thinning ', nc0, "-->", nc
!//////////////
            endif
         endif
      enddo
      end
!     *******************************
      subroutine mvdata(f, j, nc)
      implicit none
      include "Zprivate2.f"

      integer j
      real*4  f(maxp, n)

      real*4 temp(maxp)
      integer i, nc, l

      nc = 0

      do i = 1, idx(j)
         l = si(i)
         if(l .gt. 0 ) then
            nc = nc + 1
            temp(nc) = f(l, j)
         endif
      enddo

      do i = 1, nc
         f(i, j) = temp(i)
      enddo


      end
      subroutine mvdatai(f, j, nc)
      implicit none
      include "Zprivate2.f"
      integer j
      integer*2  f(maxp, n)

      integer*2  temp(maxp)
      integer i, nc, l
      nc = 0

      do i = 1, idx(j)
         l = si(i)
         if(l .gt. 0 ) then
            nc = nc + 1
            temp(nc) = f(l, j)
         endif
      enddo

      do i = 1, nc
         f(i, j) = temp(i)
      enddo


      end

      subroutine mkfiles
      implicit none
      include "Zprivate2.f"

      character* 130  filename

      integer  filec
      integer  klena
      data filec/0/
      save filec

      integer i, j

      do j = jmin, jmax
         if(idx(j) .gt. 0) then
            filec = filec + 1
            write(filename,'(a, a, i5.5,a)') dir(1:klena(dir)),
     *      "/ts", j, ".dat"
!     *      "/ts", filec, ".dat"
            open(20, file=filename, form='formatted')
            do i= 1, idx(j)
               write(20,'(3g16.8,i3,i3)') 
     *           x(i, j), y(i,j), z(i,j), codex(i,j), chgx(i,j)
            enddo
            close(20)
         endif
      enddo
      end

      subroutine mkskelfiles
      implicit none
      include "Zprivate2.f"

      character* 130  filename

      integer  filec
      integer  klena
      data filec/0/
      save filec

      integer i, j

      do j = jmin, jmax
         if(idx(j) .gt. 0) then
            filec = filec + 1
            write(filename,'(a, a, i5.5, a)') dir(1:klena(dir)),
     *      "/ts", j,  ".skel"
!     *      "/ts", filec, ".skel"
            open(20, file=filename, form='formatted')

            write(20, '(a)') "SKEL"
            call wtaskel(20, x(1,j), y(1,j), z(1,j), 
     *       codex(1,j), chgx(1, j), idx(j), mulalpha)

            close(20)
         endif
      enddo
      end

      subroutine mkafile
      implicit none
      include "Zprivate2.f"

      integer i, j
      character* 130  filename

      integer filec
      integer klena
      data filec/0/
      save filec

      if( filec .eq. 0 ) then
         filename = dir(1:klena(dir))//"/timesorted.dat"
         open(20, file=filename,
     *            form='formatted')
         filec = 1
      endif

      do j = jmin, jmax
         do i= 1, idx(j)
            write(20,'(3g16.8, i3,i3)' ) 
     *        x(i, j), y(i,j), z(i,j), codex(i,j), chgx(i,j)
         enddo
         if( idx(j) .gt. 0 ) write(20,*)
      enddo
      end
      subroutine mkaskelfile
      implicit none
      include "Zprivate2.f"

      integer i, j
      character* 130  filename

      integer filec
      integer klena
      data filec/0/
      save filec

      if( filec .eq. 0 ) then
         filename = dir(1:klena(dir))//"/timesorted.skel"
         open(20, file=filename,
     *            form='formatted')
         filec = 1
      endif

      do j = jmin, jmax
         if(idx(j) .gt. 0) then
            write(20,'(a)') "SKEL"
            call wtaskel(20,
     *      x(1, j), y(1,j), z(1,j), codex(1,j), chgx(1,j), idx(j), 
     *      mulalpha)
            write(20,*)
         endif
      enddo
      end
!        ****************** 
!           write a skel data for Geomview
      subroutine wtaskel(fno, xa, ya, za, code, chg,  np, mulalpha)
!           
      implicit none
      integer fno, np
      real*4 xa(np),   ya(np), za(np)
      real*4 mulalpha
      integer*2 code(np), chg(np) 

      real  r, g, b, alpha
      integer i
      if(np .gt. 0) then
         write(fno, '(2i9)') np, np
         do i = 1, np
            write(fno, '(3g16.8)' ) xa(i), ya(i), za(i)
         enddo
         do i= 1, np
            call code2rgb(code(i), chg(i), r, g, b, alpha)
            write(fno, '("1 ", i8, 3f5.2,f6.3)')
     *       i-1, r, g, b, alpha*mulalpha
         enddo
      endif
      end
      subroutine code2rgb(codei, chgi, r, g, b, alpha)
      implicit none
      integer*2  codei, chgi
      real r, g, b, alpha

      integer ncolor, mncolor
      parameter ( mncolor = 17 ) 
      type colortab
         integer*2 code
         integer*2 chg
         real      r
         real      g
         real      b
         real      alpha
      end type colortab
      type(colortab):: tab(mncolor)

      character*80 input
      integer i

      integer first
      data first/0/
      save first, ncolor

      
      if(first .eq. 0) then
         input = ' '
         open(13, file='colortab')
         i = 0
         do while(.true.)
           read(13,'(a)',end=10) input
           if(input(1:1) .ne. "#" .and. input .ne. ' ') then
              i = i + 1
              if(i .gt. mncolor)  then
                 write(0,*) ' too many color spec. in colortab'
                 stop 99999
              endif
              read(input,*)
     *          tab(i).code, tab(i).chg, 
     *          tab(i).r, tab(i).g, tab(i).b,
     *          tab(i).alpha
           endif
         enddo
 10      continue
         ncolor = i
         first = 1
         close(13)
      endif

      do  i = 1, ncolor
         if(tab(i).code .eq. codei .and. tab(i).chg .eq. chgi) goto 100
      enddo
      i = ncolor
 100  continue
      r = tab(i).r
      g = tab(i).g
      b = tab(i).b
      alpha = tab(i).alpha
      end

            
      
