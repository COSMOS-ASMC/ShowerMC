!         mufug: follow muon under ground (for 1ry e>500 tev, 80 hz,
!                for frejus site)
!      **************************************************************
!      *  usage:  if you want to run this program at tss terminal,
!      *          0)    use mujcl0 to make load module library in
!      *                #load.load. (once in a day).
!      *          1)    mylib lib1(#load.load)
!      *          2)    alloc f(sysinc) da('c2g5100.cosmos.gem') shr
!      *                 1) and 2) are needed for each session.
!      *          3)    change parameters such as  described later,
!      *                if needed.
!      *          4)    issue the fort command.
!      *     --------------------------------------------------
!      *          if you want to run this one at batch mode,
!      *          0)  same as previous 0)
!      *          1)  same as 3) previous
!      *          2)  use mufugjcl.  if the job dosen't finish,
!      *              you must see the final shower # and input
!      *              it to laste in the next run.
!      *
!      *  at run time, you can change the following parameters.
!      *    input:  input dataset name  (output from cosmos).
!      *    output: output dataset name (input to mudet. vbs binary)
!      *            if this is empty at the time of cont job (i.e.,
!      *            laste > 0),  program will fail. (if such strage
!      *            case is needed,
!      *            drop open(8,..) and use ft08 with
!      *            mod.)
!      *    rockds: dataset name containing rock profile
!      *    jobt:   cpu time alloted to the job (sec)
!      *    laste:  last shower # processed in the previous job
!      *    nend:   final shower # to be processed.
!      *   ----------------------------------------------------------
!      *
!      *  the following may be changed depending on the situation.
!      *  they can be changed only by modifiing the this file.
!      *  (see mufug0).
!      *  emin: minimum energy (tev) of muon to be followed in
!      *        the rock. (default is 400 mev=kitentic energy)
!      *  mode: 1 or 3 or 4 depending the 1,3,4-dimensional simulation
!      *        if only range fluctuation is needed, 1 may be used.
!      *        (default is 3).
!      *  epsp: fractional energy loss below this value is treated
!      *        as a continuous loss.  (pair creation)
!      *        (default is 1/10.)
!      *  epsb: same as epsp for brems.
!      *        (default is 1/20.)
!      *  epsh: same as epsp for hadron production.
!      *        (default is 1/20.)
!      * -----------------------------------------------------------
!      *
!      *  note for installing mufug on other machines.
!      *    following is machin dependent features.
!      * 1) call  timei and call timec.  these are subroutines to
!      *    check cpu time availability. may be dropped if the user
!      *    can use cpu time freely.
!      *
!      * module names needed for exportation.
!      *        $mucom
!      *        $mufol
!      *        $mufug
!      *        mubrem
!      *        mucset
!      *        mudedx
!      *        mufolw
!      *        nufug
!      *        muhadr
!      *        mumin
!      *        mupair
!      *        musub
!      *        rnd
!               (timei)
!      *
!      *
!      ****************************************************************
!        next external is needed when name option is used in
!       the compilation of subroutines
       external blockd
!            init
       call mufug0
!       *** until loop*** 
       do while (.true.)
          call murd1(7,  icon)
          if(icon .eq. 0) then
              call murd2(7)
!                  underground development
              call mufug
          endif
          need=1
          call timec(need, jcon)
       if         (icon .ne. 0 .or. jcon .ne. 0)
     *                    goto 1000
       enddo
 1000  continue
!          end of all events
       call mufuge(jcon)
      end
!         init mufug.
      subroutine mufug0
       include  'Zmufug.f'
         data emin/.5e-3/, mode/3/, epsp/.1/, epsb/.05/,
     *        epsh/.05/
         parameter (pi=3.141592, Torad=pi/180.)
!
         character*32 rockds/'c2g5100.rock.data(frejus)'/
         character*24 input/'c2s5042.#frj2.data'/
         character*24 output/'c2s5042.#gd2.data'/
         namelist /contn/ rockds, jobt, laste, nend, input,output
!
       laste=0
       nend=9999999
       jobt=55
       write(*,*) ' cont information '
       read(*,*) jcont
       if(jcont .eq. 1) then
            read(9, contn)
       else
           write(*,*) ' enter rockds, jobt, laste, nend (='
           write(*,*) rockds, jobt, laste, nend, ' )'
           read(*,*) rockds, jobt, laste, nend
           write(*,*) rockds, jobt, laste, nend
!
          write(*,*) ' enter input ds(=',input, ') output ds(=',output,
     *            ')'
           read(*,*) input, output
           write(*,*) input, output
       endif
!
!        original data.
       open(7, file=input, action='read',
     *         status='shr')
!        output data.
       open(8, file=output,action='both')
!
!         init timer
       call timei('00.00.00', '00.00.00', jobt)
!         open rock profile file
       call muroc0(rockds, a, z, zba, z2ba, rho, bh, bv, beta)
       cosb=cos(beta*Torad)
       sinb=sin(beta*Torad)
       call mucset(z, a, zba, z2ba)
       call mupair(epsp)
       call mubrem(epsb)
       call muhadr(epsh)
       call muradl(z, zba, z2ba,   x0ing)
       write(*,*) ' x0ing=',x0ing
       call mufol0(rho, x0ing, emin,  mode)
       if(laste .gt. 0) then
!             this should be cont job. output should be
!           appended.  skip existing data
          call muskip(8)
       endif
       return
!      *******************
       entry mufuge(jcon)
!      *******************
!          end of all events or end of time
          write(*,*) ' last shower # processed=', nshwno
          if(jcon .eq. 0) then
             write(*,*) ' all events ended '
          else
             write(*,*) ' process to be continued'
             write(*,*) ' give  1 to jcont '
             laste=nshwno
             write(9, contn)
          endif
       end
      end
      subroutine muskip(io)
           do   i=1, 9999999
               read(io, end=900)
           enddo
  900     continue
          backspace io
      end
      subroutine mufug
!         ------------------------------------------------------
       include  'Zmufug.f'
           real*8 tx, ty, tz
           data twopi/6.283184/
!             obtaine incident particle direction in detector system
           wx=w1inp*cosb + w2inp*sinb
           wy=-w1inp*sinb + w2inp*cosb
           wz= w3inp
!             get teta and fai(in rad) teta, fai in real*4
           call mudtoa(wx, wy, wz, teta, fai)
           if(fai .lt. 0.) then
               fai=fai+twopi
           endif
!             get path length in g/cm**2
           call murock(teta, fai, pathg)
!              get minimum energy to reach the depth
           call mumin(pathg, ec)
!              inform depth to muon follower
           call mufol1(pathg, 1)
!              # of muon counter to reach the depth
           muc=0
            do   n=1, nmu
               if( oba(3, n) .gt. ec) then
                    e=oba(3,n)
                    x=oba(1,n)
                    y=oba(2,n)
                    tx=aoba(1,n)
                    ty=aoba(2,n)
                    tz=sqrt(1.d0- (tx**2+ty**2))
                    call mufolw(e, x, y, tx, ty, tz,
     *              eo, xo, yo, txo, tyo, tzo, m)
                    if(m .gt. 0) then
!                        m=1 because  1 is given in mufol1
!                        save the information
                       muc=muc+1
                       eats(muc)=e
                       ed(muc)=eo
                       xd(muc)=xo
                       yd(muc)=yo
                       txd(muc)=txo
                       tyd(muc)=tyo
                       tzd(muc)=tzo
!$$$$$$$$$$$$$$$$$
!                      write(*,*) ' reach; e=',eo, ' x=',xo,
!    *                  ' y=',yo, ' tx=',txo, ' ty=',tyo,' tz=',tzo
!                   write(*,*) ' e, x, y, tx, ty, tz=', e, x, y, tx, ty,
!    *              tz
!                   else
!$$$$$$$$$$$$$$$$$
                    endif
                endif
            enddo
            if(muc .gt. 0) then
                call muout(muc)
            endif
      end
!        convert coord. in 1ry system into detector system.
      subroutine muout(muc)
       include  'Zmufug.f'
!                intit for 1ry to * system converstion matrix
           call  mugeo0(bh, bv, w1inp, w2inp, w3inp, tm)
!             convert to detector system
            do   i=1, muc
               call mucvds(xd(i), yd(i), txd(i), tyd(i), tzd(i),
     *            xd(i), yd(i), txd(i), tyd(i), tzd(i))
!$$$$$$$$$$$$$$
!        write(*,*) ' ed=',ed(i), ' xd=',xd(i), ' yd=',yd(i),
!    *        ' txd=',txd(i), ' tyd=',tyd(i), ' tzd=',tzd(i),
!    *        ' w3inp=',w3inp
!$$$$$$$$$$$$$
            enddo
           write(8) nshwno, e1ry, ihg1ry, sngl(wx), sngl(wy),
     *     sngl(wz), zfirst
           write(8)
     *     muc, (eats(i), ed(i), xd(i), yd(i), txd(i), tyd(i),
     *              tzd(i),  i=1, muc)
       end
       subroutine mucvds(xin, yin, txin, tyin, tzin, xo, yo,
     *      txo, tyo, tzo)
       include  'Zmufug.f'
!           convert x,y in 1ry system into detector system;
!           z is made to be on the same level.
!             adjust x,y so that ptcl reach the same z* plane
!             ptcl should be moved by el
!              angle conversion
           call muptos(txin, tyin, tzin, txs, tys, tzs)
           el= - (tm(3,1)*xin+ tm(3,2)*yin)/tzs
!            x,y in * system
           call muptos(xin+el*txin, yin+el*tyin, el*tzin, xs, ys, zs)
!             rotation by beta
           xo=xs*cosb + ys*sinb
           yo=-xs*sinb+ ys*cosb
           txo=txs*cosb+tys*sinb
           tyo=-txs*sinb+tys*cosb
           tzo=tzs
       end
!
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
             write(*,*) ' neglected '
             write(*,*) '*****************************'
         endif
       end
