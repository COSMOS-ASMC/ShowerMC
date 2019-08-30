#define G95
!  #define ABSOFT
      implicit none
#include "Zprivate2.f"

      character*120  tracefile, timeprofile, filename
      integer  i
      integer klena, ntime
      integer ntimem
      parameter (ntimem = 10)
      real ta1(ntimem), ta2(ntimem), dta(ntimem), 
     *  pixela(ntimem), rmaxa(ntimem),
     *  maxppta(ntimem), maxthina(ntimem), 
     *  offseta(ntimem)
      character*120  dira(ntimem)

      do i = 1, 9
         codesel(i) = 0
      enddo
      read(*,*) 
      read(*,*)  
     *  tracefile, timeprofile, split, outtype
      read(*,*)
      read(*,*)  chgsel,  codesel
!
!       read det to primary system conversion matrix
      read(*,*)
      read(*,*) aax, aay, aaz
      read(*,*) bbx, bby, bbz
      read(*,*) ccx, ccy, ccz
                  

      ncodes = 0
      do i = 1, 9
         if( codesel(i) .eq. 0) goto 10
         ncodes = ncodes + 1
      enddo
 10   continue
      select = 0
      if(chgsel .ne. -1)  select = select + 1
!         code selection imposed
      if(ncodes .gt. 0)  select = select + 2
!
!
      write(0,*)  tracefile(1:klena(tracefile)), " ",
     *  timeprofile(1:klena(timeprofile)), split, outtype
      write(0,*)  chgsel, codesel
      write(0,*) aax, aay, aaz
      write(0,*) bbx, bby, bbz
      write(0,*) ccx, ccy, ccz
                  
#ifndef ABSOFT
      call flush(0)
#endif
      open(11, file=timeprofile(1:klena(timeprofile)) )
      ntime = 0
      write(0,*) ' profile opened'
      read(11,*) 
      do while(.true.)
         i = ntime +1
         read(11,*, end=200) ta1(i), ta2(i), dta(i), pixela(i), 
     *   rmaxa(i),  maxppta(i), maxthina(i), 
     *   offseta(i),  dira(i)
         if( ta1(i) .eq. 0. .and. ta2(i) .eq. 0.) goto 200

         write(0, *) ta1(i), ta2(i), dta(i), pixela(i), 
     *   rmaxa(i), maxppta(i), maxthina(i), 
     *   offseta(i),  dira(i)(1:klena(dira(i)))

         ntime = ntime +1
         if(rmaxa(i) .eq. 0.) then
            rmaxa(i) = 1.e10
         endif
      enddo
 200  continue
      close(11)
      do i = 1, ntime
         filename = ' '
         filename = dira(i)(1:klena(dira(i)))//"/dummy"
         open(11, file=filename,  err=210) 
         close(11)
      enddo
      goto 220
 210  continue
      write(0,*) ' following dir may be missing:= ',
     *   dira(i)(1:klena(dira(i)))
      stop 3333


 220  continue
      open(11, file=tracefile(1:klena(tracefile)), 
     *         form='formatted', err=230 ) 
      goto 240
 230  continue
      write(0,*) ' trace file may be missing:= ',
     *  tracefile(1:klena(tracefile))
      stop 4444
 240  continue
      write(0,*) ' trace file opend'
      do i = 1, ntime
         if(maxppta(i) .eq. 0) then
            maxppt = maxp
         else
            maxppt = maxppta(i) 
         endif
         maxthin = maxthina(i)
         offset = offseta(i)
         dir = dira(i)
         call slice1tbin(
     *    ta1(i), ta2(i), dta(i), pixela(i), rmaxa(i) )
         rewind 11
      enddo
      end

      subroutine slice1tbin(mint, maxt, dt, pixelin, rmax)
      implicit none
#include "Zprivate2.f"
      real*4  mint, maxt, dt, t1, t2, rmax, pixelin

      rmax2 = rmax**2
      t2 = mint
      pixel = pixelin
      do while (t2 .lt. maxt)
         t1 = t2
         t2 = min(t1+ dt*n , maxt)
         jmax = 0
         jmin = 10000000

         call slice(t1, t2, dt)
         write(0,*) ' t1=',t1, 't2=',t2, ' slice ended'
         write(0,*) ' maxppt=', maxppt
         call thinning(jmin, jmax)
         write(0,*) ' thinning ended'
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
#include "Zprivate2.f"

      real*4  mint, maxt, dt
      integer i, j, kk
      real*4  t,  k
      integer*2 code, chg 
      integer it1, it2, kkk
      real*4  xx1, yy1, zz1, time1
      real*4  xx2, yy2, zz2, time2
      real*4  xxp, yyp
      real*4  kerg1, kerg2
      logical chkcode
      character* 120 input
      real*8 linec
      data linec/0.0/
      save linec
      
      chkcode = ncodes .gt. 0

      do j = 1, n
         idx(j) = 0
         thinc(j) = 0
      enddo
      
      do while (.true.)

         input = ' '
         read(11, '(a)', end=1000 )  input
         linec = linec + 1.d0
         if( mod(linec, 1.d7) .eq.  0 ) then
            write(0,*)
     *      ' ------------------------------------linec=',linec
#ifndef ABSOFT
            call flush(0)
#endif
         endif
         if( input .ne. " " ) then
            read(input, * ) xx1, yy1, zz1,   code, kerg1, chg,
     *         time1
            if(select .eq. 0) goto 100
!                  select all
            if( chgsel .ne. -1 )  then
!                  charge selection
               if(chgsel .eq. 0 .and. chg .ne. 0) goto 950
               if(chgsel .eq. 1 .and. chg .eq. 0) goto 950
            endif
            if( chkcode ) then
               do kkk = 1, ncodes
                  if(codesel(kkk) .eq. code) goto 100
               enddo
               goto 950
            endif
!              
 100        continue
            do while( input .ne. " " ) 
               read(11, '(a)', end =1000)  input
               linec = linec + 1.d0
               if( mod(linec, 1.d7) .eq.  0 ) then
                  write(0,*)
     *            ' ------------------------------------linec=',
     *            linec
#ifndef ABSOFT
                  call flush(0)
#endif
               endif

               if( input .ne. " " ) then
                  read(input, * ) xx2, yy2, zz2,  code, kerg2, chg,
     *            time2 
                  xxp = xx2*aax + yy2*aay   ! aaz =0
                  yyp = xx2*bbx + yy2*bby + zz2*bbz
                  it1 = (time1- mint)/dt +1
                  it2 = (time2- mint)/dt +1

                  if(xxp**2 + yyp**2 .gt. rmax2) goto 900
                  if(it1 .lt. 1) goto 900
                  if(it2 .gt. n) goto 900
                  if(time2 .gt. maxt) goto 900
                  jmin = min(jmin, it1)
                  jmax = max(jmax, it2)
                  do  j = it1, it2
                     t = j*dt + mint
                     if( idx(j) .eq. maxppt  .and. 
     *                  thinc(j) .lt. maxthin ) then
                        write(0,*)
     *                    'ptcls at time=', t , ' > ',  maxppt
                        write(0,*)
#ifndef ABSOFT
                        call flush(0)
#endif
                        call thinning(j, j)
                        thinc(j) = thinc(j) + 1
                     endif
                     idx(j) = idx(j) + 1
                     if(idx(j) .le. maxppt) then
!                           (x2-x1)/(t2-t1) *(t-t1) +x1
                        if(time2 .eq. time1) then
                           k =0.
                        else
                           k = (t - time1) /(time2-time1)
                        endif
                        x(idx(j), j) = (xx2-xx1) * k + xx1
                        y(idx(j), j) = (yy2-yy1) * k + yy1
                        z(idx(j), j) = (zz2-zz1) * k + zz1
                        codex(idx(j), j) = code
                        chgx(idx(j), j) = chg
                     endif
                  enddo   
 900              continue
                  xx1 = xx2
                  yy1 = yy2
                  zz1 = zz2
                  time1 = time2
               endif
            enddo
         endif
 950     continue
      enddo
 1000 continue
      end

!     ***********************************
!      do thinning if two or more data fall in the dame pixel.
!     **********************************
      subroutine thinning(i1, i2)
      implicit none
#include "Zprivate2.f"
!        culling overlapping points in a pixel
!           this is very rough thinning
!     sorr each x
      integer i1, i2
      integer nc, j, i, k1, k2, l, nc0
      logical  dothin
      integer ntp
      

      do j = i1, i2
         ntp= min(idx(j), maxppt)
         dothin = .false.
         if( ntp  .gt. maxppt/2.0 ) then
            call kqsortr(x(1,j), si, ntp)
            do i = 1, ntp-1
               k1 = si(i)
               if( k1  .lt. 0 ) goto 20
               do l = i+1,  ntp
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
               nc0 = ntp
               call mvdata(x, j, nc)
               call mvdata(y, j, nc)
               call mvdata(z, j, nc)
               call mvdatai(codex, j, nc)
               call mvdatai(chgx, j, nc)  
               write(0,*)  ' at time idx=',j,
     *          ' thinning ', nc0, "-->", nc
#ifndef  ABSOFT
               call flush(0)
#endif
               idx(j) = nc
            endif
         endif
      enddo
      end
!     *******************************
      subroutine mvdata(f, j, nc)
      implicit none
#include "Zprivate2.f"

      integer j
      real*4  f(maxp, n)

      real*4 temp(maxp)
      integer i, nc, l, ntp


      nc = 0
      ntp = min(idx(j), maxppt)
      do i = 1,  ntp
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
#include "Zprivate2.f"
      integer j
      integer*2  f(maxp, n)

      integer*2  temp(maxp)
      integer i, nc, l, ntp


      nc = 0
      ntp = min(idx(j), maxppt)

      do i = 1, ntp
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
#include "Zprivate2.f"

      character* 130  filename
      integer filec
      integer  klena
      integer i, j
      data filec/0/
      save filec

      do j = jmin, jmax
         if(idx(j) .gt. 0) then
            filec = filec +1
            write(filename,'(a, a, i5.5,a)') dir(1:klena(dir)),
     *      "/ts", filec+offset, ".dat"
            open(20, file=filename(1:klena(filename)), 
     *        form='formatted')
            do i= 1, min( idx(j), maxppt)
               write(20,'(3g16.8,i3,i3)') 
     *           x(i, j), y(i,j), z(i,j), codex(i,j), chgx(i,j)
            enddo
            close(20)
         endif
      enddo
      end

      subroutine mkskelfiles
      implicit none
#include "Zprivate2.f"

      character* 130  filename
      integer filec

      integer  klena, ntp

      integer i, j
      data filec/0/
      save filec


      write(0,*) ' jmin,max=', jmin, jmax


      do j = jmin, jmax
         ntp = min(idx(j), maxppt)
         if(idx(j) .gt. 0) then
            filec = filec + 1
            write(filename,'(a, a, i5.5, a)') dir(1:klena(dir)),
     *      "/ts", filec+offset,  ".skel"

            open(20, file=filename(1:klena(filename)),
     *         form='formatted')
            write(20, '(a)') "SKEL"
            call wtaskel(20, x(1,j), y(1,j), z(1,j), 
     *       codex(1,j), chgx(1, j), ntp)

            close(20)
         endif
      enddo
      end

      subroutine mkafile
      implicit none
#include "Zprivate2.f"

      integer i, j
      character* 130  filename

      integer filec
      integer klena
      data filec/0/
      save filec

      if( filec .eq. 0 ) then
         filename = dir(1:klena(dir))//"/timesorted.dat"
         open(20, file=filename(1:klena(filename)),
     *            form='formatted')
         filec = 1
      endif

      do j = jmin, jmax
         do i= 1, min(idx(j), maxppt)
            write(20,'(3g16.8, i3,i3)' ) 
     *        x(i, j), y(i,j), z(i,j), codex(i,j), chgx(i,j)
         enddo
         if( idx(j) .gt. 0 ) write(20,*)
      enddo
      end
      subroutine mkaskelfile
      implicit none
#include "Zprivate2.f"

      integer i, j
      character* 130  filename

      integer filec
      integer klena, ntp
      data filec/0/
      save filec

      if( filec .eq. 0 ) then
         filename = dir(1:klena(dir))//"/timesorted.skel"
         open(20, file=filename(1:klena(filename)),
     *            form='formatted')
         filec = 1
      endif

      do j = jmin, jmax
         ntp = min(idx(j), maxppt)
         if(idx(j) .gt. 0) then
            write(20,'(a)') "SKEL"
            call wtaskel(20,
     *      x(1, j), y(1,j), z(1,j), codex(1,j), chgx(1,j), ntp)
            write(20,*)
         endif
      enddo
      end
!        ****************** 
!           write a skel data for Geomview
      subroutine wtaskel(fno, xa, ya, za, code, chg,  np)
!           
      implicit none
      integer fno, np
      real*4 xa(np),   ya(np), za(np)
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
            write(fno, '("1 ", i8, 3f5.2)')
     *       i-1, r, g, b
         enddo
      endif
      end
      subroutine code2rgb(codei, chgi, r, g, b, alpha)
      implicit none
      integer*2  codei, chgi
      real r, g, b, alpha

      integer ncolor, mncolor
      parameter ( mncolor = 17 ) 
!      structure /colortab/
         integer*2 tcode(mncolor)
         integer*2 tchg(mncolor)
         real      tr(mncolor)
         real      tg(mncolor)
         real      tb(mncolor)
         real      talpha(mncolor)
!      end structure
!      record /colortab/ tab(mncolor)

      character*80 input
      integer i

      integer first
      data first/0/
      save first, ncolor, tcode, tchg, tr, tg, tb, talpha


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
     *          tcode(i), tchg(i), 
     *          tr(i), tg(i), tb(i),
     *          talpha(i)
           endif
         enddo
 10      continue
         ncolor = i
         first = 1
         close(13)
      endif

      do  i = 1, ncolor
         if(tcode(i) .eq. codei .and. tchg(i) .eq. chgi) goto 100
      enddo
      i = ncolor
 100  continue

      r = tr(i)
      g = tg(i)
      b = tb(i)
      alpha = talpha(i)
      end

            
      
