!         musurf: analyse muons at surface of the earth.
!      #gd.data: fmt=3, 111 # of data =3  (use nmu command)
!      #gd2.data: fmt=3, 111 # of data=4  (use muer command)
!      **************************************************************
!      *  usage:  if you want to run this program at tss terminal,
!      *          1)    alloc f(sysinc) da('c2g5100.cosmos.gem') shr
!      *          2)    change parameters such as  described later,
!      *                if needed.
!      *          3)    issue the fort command.
!      *     --------------------------------------------------
!      *
!      *  at run time, you can change the following parameters.
!      *    input:  input dataset name  (output from cosmos).
!      *    outpt1: output dataset name (input to gd)
!      *            if this is empty at the time of cont job (i.e.,
!      *            laste > 0),  program will fail.
!      *       for quantities: 1 data /event
!      *    outpt2: output dataset name (input to gd)
!      *            if this is empty at the time of cont job (i.e.,
!      *            laste > 0),  program will fail.
!      *       for quantities:  nmu data/event
!      *    jobt:   cpu time alloted to the job (sec)
!      *    laste:  last shower # processed in the previous job
!      *    nend:   final shower # to be processed.
!      *   ----------------------------------------------------------
!      *
!      *  note for installing musurf on other machines.
!      *    following is machin dependent features.
!      * 1) call  rndc(u)  which should give random # in 0 to 1.0
!      * 2) call  timei and call timec.  these are subroutines to
!      *    check cpu time availability. may be dropped if the user
!      *    can use cpu time freely.
!      *
!      *
!      ****************************************************************
       logical first
        do   kkk=1, 100
          first=.true.
!            init
          call musur0
!          *** until loop*** 
          do while (.true.)
             call murd1(7,  icon)
             if(first) then
                 first=.false.
                 call mucap
             endif
             if(icon .eq. 0) then
                 call murd2(7)
!                    analysis at surface
                 call musurf
             endif
             need=1
             call timec(need, jcon)
          if         (icon .ne. 0 .or. jcon .ne. 0)
     *                       goto 1000
          enddo
 1000     continue
!            end of all events
          call musure(jcon)
        enddo
      end
!         init musurf
      subroutine musur0
       include  'Zmufug.f'
       common /$/ emin
       common /$$/ intm
       character*4 intm
!
         character*24 input/'c2s5001.#data.data'/
         character*24 outpt1/'c2s5001.#gd.data'/
         character*24 outpt2/'c2s5001.#gd2.data'/
         character*16 cap1,cap2, cap3, cap4, cap5, cap6
         character*70 ttl
         data
     *   cap1/'r(from axis)'/,cap2/'energy(tev)'/,cap3/'1st z(g/cm2)'/
         data cap4/'# of muons'/, cap5/'r(from g.c)'/
         data cap6/'# of z orig'/
!
       laste=0
       nend=9999999
       jobt=55
       emin=1.
       write(*,*) ' enter emin, jobt, laste, nend (='
       write(*,*) emin, jobt, laste, nend, ' )'
       read(*,*)  emin, jobt,  laste, nend
       write(*,*) emin,  jobt, laste, nend
       if(emin .eq. 0.) then
           stop
       elseif(emin .eq. 1.) then
           intm='lund'
           modd=0
       elseif(emin .ne. 1.) then
           modd=1
       endif
       if(modd .eq. 0) then
!
!        write(*,*) ' enter input ds(=',input, ') outpt1 ds(=',outpt1,
!    *   ') output2 ds(=',outpt2,')'
!        read(*,*) input, outpt1, outpt2
         write(*,*) input, outpt1, outpt2
         write(*,*) ' enter int model(=', intm, ')'
         read(*,'(a)') intm
         if(intm .eq. '/' .or. intm .eq. '   ') then
             intm='lund'
         endif
       endif
!
!        original data.
       open(7, file=input, action='read',
     *         status='shr')
!        output data.
       open(8, file=outpt1,action='both')
       open(9, file=outpt2,action='both')
!
!         init timer
       call timei('00.00.00', '00.00.00', jobt)
       if(laste .gt. 0 .or. modd .gt. 0 ) then
!             this should be cont job. output should be
!           appended.  skip existing data
          call muskip(8)
          call muskip(9)
       endif
       return
!*******************
      entry mucap
!******************
       write(ttl,
     * '('' intm='',a,''e0='',f7.0,''tev:'',a,'' cos='',f4.2,
     * '' emin='',f4.1,'' tev'')') intm, e1ry, k1ry, w3inp, emin
       if(laste .eq. 0) then
             write(8) ttl
             write(9) ttl
             write(8) cap3, cap4, cap6
             write(9) cap2, cap4, cap5, cap1
       endif
      end
      subroutine muskip(io)
           do   i=1, 9999999
               read(io, end=900)
           enddo
  900     continue
          backspace io
      end
      subroutine musurf
!         ------------------------------------------------------
       include  'Zmufug.f'
       common /$/ emin
       common /$$/ intm
       character*4 intm
           dimension rcnt(800)
!          wx=w1inp
!          wy=w2inp
!          wz=w3inp
           nmux=0
            do   i=1, nmu
               if(oba(3, i) .ge. emin) then
                   nmux=nmux+1
                   call kmover(oba(1,i), 1, 10, oba(1, nmux), 1)
               endif
            enddo
           nmu=nmux
           call murcnt(oba(1,1), oba(2,1), 10, nmu, rcnt)
            do   n=1, nmu
                    e=oba(3,n)
                    x=oba(1,n)
                    y=oba(2,n)
!                   znnc=oba(4,n)
!                   znnc=oba(5,n)
                    rax=sqrt(x**2+y**2)
                write(9)  e, float(nmu), rcnt(n), rax
            enddo
           call mufdz(oba(4,1), 10, nmu, ndifz)
           write(8)  zfirst, float(nmu), float(ndifz)
      end
!     ============================================================
       subroutine murd1(mtno,  icon)
!         icon=0   next is data
!              1   no more data
!
       include  'Zmufug.f'
!
          dimension ir1(2),  seobp(9, 1), nobtp(9, 1),
     *           anobp(9, 1), rsvd1(10)
!    *           anobp(9, 1), rsvd1(10), suftbl(16)
          data j2/1/, kindma/9/, mnlvl/1/, l/1/
!         used to identify the top of 1 shower on mt.
          data sttmk/z7ffffff1/
!           //             top of 1 level of ec1 type
          data lvlmk/z7ffffff2/
          icon=0
!          *** until loop*** 
          do while (.true.)
!              *** until loop*** 
              do while (.true.)
                  read(mtno, end=8900) topcd
              if         (topcd .eq. sttmk)
     *                           goto 510
              enddo
  510         continue
!
              read(mtno, end=8900) nshwno, ir1, ir2, e1ry, k1ry,
     *        ( (seobp(i,j),i=1, kindma),j=1,j2),
     *        ( (nobtp(i,j),i=1, kindma),j=1,j2),
     *        ( (anobp(i,j),i=1, kindma),j=1,j2),
     *        w1inp, w2inp, w3inp, ihg1ry,
     *        rsvd1, nthnc, zfirst
!    *        ,suftbl
!    *        ,(znthnc(i),anthnc(i),enthnc(i),nnpa(i),
!    *        ettma(i),i=1,nthnc)
          if         (nshwno .gt. laste)
     *                       goto 555
          enddo
  555     continue
          if(nshwno .gt. nend) then
            goto 8900
          endif
!
      loop=0
      ll=0
       do   while (loop .le. mnlvl .and. ll .ne. l)
!         *** until loop*** 
         do while (.true.)
            read(mtno, end=8900)lvlmkx
         if         (lvlmkx .eq. lvlmk)
     *                      goto 210
         enddo
  210    continue
!
         read(mtno, end=8900) ll, ntotal, wtotal, setotl
         loop=loop+1
       enddo
      if( ll  .ne. l) then
         write(*,*)' level=',l,' not exists'
         stop
      endif
      return
!
 8900 continue
      icon=1
      return
      end
      subroutine murd2(mtno)
       include  'Zmufug.f'
!          k= 4---> mu,  8--> neu e   9---> neu mu
!          np:  # of ptcl
         parameter (muon=4)
!         *** until loop*** 
         do while (.true.)
            read(mtno) nk, wno, k,l,np, ncum, npg, ntpg,
     *                   ((oba(i,j),i=1,10),   j=1, np)
!ccc        if(obangl) then
                 read(mtno) np, ((aoba(i,j),i=1,2),j=1,np)
!ccc        endif
         if         (k .eq. muon)
     *                      goto 100
         enddo
  100    continue
         nmu=np
         if(ntpg .gt. 1) then
             read(mtno) nk, wno, k, l, np, ncum, npg, ntpg,
     *           ((oba(i,j), i=1, 10), j=nmu+1, nmu+np)
             read(mtno) np, ((aoba(i,j), i=1,2), j=nmu+1, nmu+np)
             nmu=nmu+np
         endif
         if(ntpg .gt. 2) then
             write(*,*) '*****************************'
             write(*,*) ' warning: > 800 muons exist but '
             write(*,*) ' muons >800 neglected '
             write(*,*) '*****************************'
         endif
       end
       subroutine musure(jcon)
       include  'Zmufug.f'
!          end of all events or end of time
          write(*,*) ' last shower # processed=', nshwno
          if(jcon .eq. 0) then
             write(*,*) ' all events ended '
             write(8) 1.e50, 1.e50, 1.e50
             write(9) 1.e50, 1.e50, 1.e50, 1.e50
          else
             write(*,*) ' process to be continued'
             write(*,*) ' give ',nshwno, ' to laste'
          endif
          close(7)
          close(8)
          close(9)
       end
       subroutine mudeco(x, y, intv, n, r, m)
          dimension x(intv,n), y(intv,n), r(*)
          if(n .gt. 10 ) then
              write(*,*) ' too many data'
          else
              m=0
               do   i=1, n-1
                   do   j=i+1, n
                      rx=sqrt(  (x(1,j)-x(1,i))**2 +
     *                (y(1,j)-y(1, i))**2 )
                      m=m+1
                      r(m)=rx
                   enddo
               enddo
          endif
       end
       subroutine murcnt(x, y, intv, n, r)
          dimension x(intv,n), y(intv,n), r(n)
           if(n .gt. 1) then
              sumx=0.
              sumy=0.
               do   i=1, n
                  sumx=sumx+x(1,i)
                  sumy=sumy+y(1,i)
               enddo
              rx=sumx/n
              ry=sumy/n
               do   i=1, n
                  r(i)=sqrt( (x(1,i)-rx)**2 +
     *                (y(1,i)-ry)**2 )
               enddo
          else
              r(1)=0.
          endif
       end
       subroutine mufdz(z, intv, n, ndif)
           dimension z(intv, n)
           if(n .le. 1) then
                ndif=n
           else
               call ssort1(z, 1, n, 4*intv,1, 4, 'r', 'a', 'all')
               ndif=1
               zx=z(1, 1)
                do   i=2, n
                   if(z(1,i) .ne. zx) then
                       ndif=ndif+1
                       zx=z(1,i)
                   endif
                enddo
           endif
       end
