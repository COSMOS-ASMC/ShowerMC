      implicit none
!   Suppose  we have n hybrid data
!     h1, h2, h3, ....hn
!   1) cp h1 h0
!   2) h0 + h2 --> h; mv  h h0
!   3) h0 + h3 --> h; mv  h h0
!   ..
!   n) h0 + hn --> h; mv  h h0
! 
!   This program add two hybrid data; h0 + hx--> h
!            environmental variable
!   file h0: HYBFILE0
!   file hx: HYBFILEX
!   file h:  HYBFILET
!

      integer ndepth
      parameter (ndepth= 50)
      real*8  ASdep(ndepth),  muunit(ndepth), sumEsize, sumEsize2
      real*8  Esize0(ndepth),
     *     age0(ndepth),  cogdep0(ndepth),
     *     SEloss0(ndepth),
     *     Ng0(ndepth), Ne0(ndepth), Nmu0(ndepth), Nhad0(ndepth),
     *     cog0, cog20
      real*8  Esizex(ndepth),
     *     agex(ndepth),  cogdepx(ndepth),
     *     SElossx(ndepth), 
     *     Ngx(ndepth), Nex(ndepth), Nmux(ndepth), Nhadx(ndepth),
     *     cogx, cog2x

      real*8  Esizet(ndepth),
     *     aget(ndepth),  cogdept(ndepth),
     *     SElosst(ndepth),
     *     Ngt(ndepth), Net(ndepth), Nmut(ndepth), Nhadt(ndepth),
     *     cogt, cog2t
            
      integer EvNo0, EvNox
      integer fn0, fnx, fnt
      integer kgetenv2
      integer klena
      integer code, subcode, charge
      integer leng, i, j, lengflesh
      integer icon0, iconx, icont
      real dd, firstz, w1, w2, w3, E0
      character*128 hyb0, hybx, hybt, flesher
      character*128 input0, inputx, inputt
      fn0 = 2
      fnx = 3
      fnt = 4
      leng = kgetenv2("HYBFILE0", hyb0)
      call copenfw2(fn0, hyb0, 1, icon0)
      if(icon0 .ne. 1)  then
         write(0,*) hyb0(1:leng)
         if( icon0 .eq. 0) then
            write(0,*) 'not exists'
         else
            write(0,*) ' cannot be opened '
         endif
         write(0,*) ' icon=',icon0
         stop 9999
      else
         write(0,*)  hyb0(1:leng), ' opened'
      endif
      leng = kgetenv2("HYBFILEX", hybx)
      call copenfw2(fnx, hybx, 1, iconx)
      if(iconx .ne. 1)  then
         write(0,*) hybx(1:leng)
         if( iconx .eq. 0) then
            write(0,*) 'not exists'
         else
            write(0,*) ' cannot be opened '
         endif
         write(0,*) ' icon=',iconx
         stop 9999
      else
         write(0,*)  hybx(1:leng), ' opened'
      endif

      
      leng = kgetenv2("HYBFILET", hybt)
      call copenfw2(fnt, hybt, 1, icont)
      if(icont .ne. 0)  then
         write(0,*) hybt(1:leng)
         write(0,*) ' cannot be opened '
         write(0,*) ' icon=',icont
         stop 9999
      else
         write(0,*)  hybt(1:leng), ' opened'
      endif
      lengflesh = kgetenv2("FLESHDIR", flesher)

      i =0
      do while(.true.)
         input0 = ' '
         read( fn0, '(a)',  end=1000 ) input0
         inputx = ' '
         read( fnx, '(a)'  ) inputx 
         if(input0 .ne. " ") then
            if(flesher(1:lengflesh) .eq. "FleshHist") then
               if(input0(1:1) .eq. "h" ) then
                  read(input0(3:klena(input0)),* ) EvNo0, code,
     *            subcode, charge, E0, w1, w2, w3,firstz, 
     *            cog0, cog20

                  read(inputx(3:klena(inputx)),*) EvNox, code,
     *            subcode, charge, E0, w1, w2, w3,firstz,
     *            cogx, cog2x

               else
                  read(input0(3:klena(input0)), *) 
     *            i, ASDep(i),  muunit(i), age0(i), cogdep0(i),
     *            Ng0(i), Ne0(i), Nmu0(i), Nhad0(i), 
     *            Esize0(i),  SEloss0(i)

                  read(inputx(3:klena(inputx)), *) 
     *            i, ASDep(i),  muunit(i), agex(i), cogdepx(i),
     *            Ngx(i), Nex(i), Nmux(i), Nhadx(i), 
     *            Esizex(i),  SElossx(i)
               endif
            else
               read(input0(3:klena(input0)), *) 
     *            i, ASDep(i),  muunit(i), age0(i), cogdep0(i),
     *            Ng0(i), Ne0(i), Nmu0(i), Nhad0(i), 
     *            Esize0(i)
               SElossx(i) =0.

               read(inputx(3:klena(inputx)), *) 
     *            i, ASDep(i),  muunit(i), agex(i), cogdepx(i),
     *            Ngx(i), Nex(i), Nmux(i), Nhadx(i), 
     *            Esizex(i)
               SElossx(i) =0.
            endif
         else
!              1 event read
            if(inputx .ne. " ") then
               write(0,*) ' event differ', inputx
               stop 1111
            endif
            do j = 1, i
               Esizet(j) = Esize0(j) + Esizex(j)
               if(Esizet(j)  .gt. 0.) then
                  aget(j) = ( Esize0(j)*age0(j) + Esizex(j)*agex(j))/
     *               Esizet(j) 
               else
                  aget(j)= 0.
               endif
               SElosst(j) = SEloss0(j) + SElossx(j)
               Ngt(j) = Ng0(j) + Ngx(j)
               Net(j) = Ne0(j) + Nex(j)
               Nmut(j) = Nmu0(j) + Nmux(j)
               Nhadt(j) = Nhad0(j) + Nhadx(j)
            enddo
            cogt =0.
            cog2t =0.
            sumEsize = 0. 
            sumEsize2 = 0. 
            do j = 1, i
               if(j .gt. 1 .and. j  .lt. i ) then
                  dd =( ASDep(j+1) - ASDep(j-1))/2.0
               elseif(j .eq. 1) then
                  dd =(ASDep(2) - ASDep(1))
               else
                  dd =(ASDep(i) - ASDep(i-1))
               endif
               cogt = cogt + Esizet(j)*ASDep(j)*dd
               sumEsize = sumEsize + Esizet(j)*dd
               
               if( aget(j) .gt.  2.0-aget(i)  )  then
                  if(j .gt. 1 .and. j  .lt. i ) then
                     dd =( ASDep(j+1) - ASDep(j-1))/2.0
                  elseif(j .eq. 1) then
                     dd =(ASDep(2) - ASDep(1))
                  else
                     dd =(ASDep(i) - ASDep(i-1))
                  endif
                  cog2t = cog2t + Esizet(j)*ASDep(j)*dd
                  sumEsize2= sumEsize2 + Esizet(j)*dd
               endif
            enddo
            cogt =  cogt / sumEsize
            if(sumEsize2 .gt. 0.) then
               cog2t =  cog2t / sumEsize2
            else
               cog2t = ASDep(i)
            endif
            do j = 1, i
               cogdept(j) = ASdep(j)/cog2t
            enddo
            write(fnt, 
     *         '("h ", i4,  3i3, 1pE11.3, 0p 3f11.7, f7.2, 2f7.0)')
     *            EvNo0, code,
     *            subcode, charge, E0, w1, w2, w3,firstz, 
     *            cogt, cog2t
            do j = 1, i
               if(flesher(1:lengflesh) .eq. "FleshHist") then
                  write(fnt, '("t ", i3, 2f7.1,  2f6.3,
     *                 1p6E11.3)')
     *                 j,
     *                 ASDep(j), muunit(j),
     *                 aget(j),  cogdept(j),
     *                 Ngt(j), Net(j), Nmut(j), Nhadt(j),
     *                 Esizet(j), SElosst(j)
               else
                  write(fnt, '("t ", i3, 2f7.1,  2f6.3,
     *                 1p5E11.3)')
     *                 j,
     *                 ASDep(j), muunit(j),
     *                 aget(j),  cogdept(j),
     *                 Ngt(j), Net(j), Nmut(j), Nhadt(j),
     *                 Esizet(j)
               endif
            enddo
            write(fnt, *) 
         endif
      enddo
 1000 continue
      write(0,*) 'all events processed '
      end


            
