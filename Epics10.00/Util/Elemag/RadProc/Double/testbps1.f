!     ****************************************************************
!     *                                                              *
!     *  testing sampling of brems and pair crossection with landau  *
!     *          and/or various corrections                          *
!     *          also compt                                          *
!     *                                                              *
!     ****************************************************************
!
!       prepare:   brems, pair and compton sampling routines
!                  give data through gd
!                  make dispav=f in histgd by permanet change
!
!
      real fmt(20), tcb(8), tcp(8), tcc(4)
      real pnt(8)/' ', ' ', ' ', '*', '+', ' ', ' ', 'x'/
!
!
!       -inc $elmag,seq=100:1800
       -inc $elmag,seq=100:1800
!
!
!
      dimension iva2(160,8),iva3(150,8), iva4(70,4)
      dimension ivo(120, 8), vo(120, 8)
      logical retry , gdmode, brem,pair
      character*4 matter/' '/
      character*100 bremc, pairc, compc
      data ecut/.1e-3/
!
      gdmode=.false.
      usees=.false.
!
!
      emass=0.511e-3
      scorec=.true.
      lcorec=.true.
      lpmx=4
      call ggdin(0)
      call txtgin
      call txtgdc(1,1,1,5,84)
      call txtglc(1,'put',1)
      call rdgin(1)
      gdmode=.true.
!
    5 continue
      call ggdclr
      call txtglc(1,'put',1)
      call txtgd(1,'give nsampl=0 to stop!',0)
!
      nsampl=10000
      e1st=10.e-3
      ir=1234567
      brem=.true.
      pair=.true.
      novlp=3
      estep=10.
!
      call txtgd(1,                             'enter 1--matter(= );
     * 2--nsampl(=10000);  3--e1st(=10.e-3);  4--ir(=1234567); 5--brem(
     *t);  6--pair(t);  7--novlp(=3);  8--estep(=10); 9--scorec(l);
     *10--lcorec(l)!', 0)
      call rdg(
     * ' !',-1,
     *                  &10, &20)
      goto 5
!
   10 continue
      read(31,*,end=15) matter, nsampl, e1st,ir,brem,pair, novlp,
     *                  estep, scorec, lcorec
   15 continue
      call elmgmd(matter, 0)
      call rdgcer(&5)
   20 continue
      if(nsampl .eq. 0) call ggdtm
      if(nsampl .eq. 0) stop
      if(.not. brem  .and. .not. pair) goto 5
      call rndc(u)
!
      novlp=min0(8,novlp)
      novlp=max0(novlp,1)
      ee=e1st
      call setsv(iva2,1,160*8,0)
      call setsv(iva3,1,150*8,0)
      call setsv(ivo,1,100*8,0)
!
!
       do   ie1=1,novlp
          if(pair)  then
               call pairt(ee,t)
               tcp(ie1)=tprob
          endif
          if(brem)  then
               call bremst(ee,ecut,t)
               tcb(ie1)=tprob
          endif
!
           do   n=1,nsampl
              if(pair) then
                  call pairt(ee,t)
                  call paire(ee,e1)
                  v=e1/ee
                  call hist(v, .5, 0.01, iva3(1,ie1),150)
              endif
              if(brem) then
                  call bremst(ee,ecut,t)
                  call bremse(ee,ecut,eg)
                  v=eg/ee
                  call histl(v, -7., .05, iva2(1,ie1), 160, vlog, jcon)
                  call hist(v, 0., .02, ivo(1,ie1), 120)
              endif
           enddo
          ee=ee*estep
       enddo
!
      if(brem) then
         bremc=' '
         write(bremc,'('' total x-section of brems'',8g9.3)')
     *   (tcb(i), i=1,novlp)
         bremc=bremc(1:26+novlp*9)//'!'
         write(6,'(1h ,a100)') bremc
      endif
      if(pair) then
         pairc=' '
         write(pairc,'('' total x-section of pair'',8g9.3)')
     *   (tcp(i), i=1,novlp)
         pairc=pairc(1:25+novlp*9)//'!'
         write(6,'(1h ,a100)') pairc
      endif
!
!
      if(brem) then
           do   i=1,novlp
!                          vf(v)=t*dn/n/dlog10v/2.3026
              call histm(iva2(1,i), 160, tcb(i)*.4343/nsampl)
              call smmul2(tcb(i)/nsampl, ivo(1,i), 1, 116,
     *        vo(1,i), 1, 'ir')
           enddo
          ee=e1st
          encode(80, 935, fmt) matter, e1st
  935     format(' matter=',a4,' fbrem(v)*v at e=',g9.3,'gev!')
!
          call histgc(bremc, .true., ' ', .false., ' ', .false.)
          call histgd(fmt,iva2(1,1), -160, -7., 0.05, 1,
     *     ' log10(v)!', 6, 6, 'dx  ', nsampl, -2)
          ee=ee*estep
          infm=-3
           do   i=2,novlp
             encode(80,1315,fmt) ee
 1315        format(1h; g9.3,' gev!')
             if(i .eq. novlp) infm=-2
             call histg2(fmt,iva2(1,i),-160,-7.,0.05,nsampl,
     *       pnt(i),infm)
             ee=ee*estep
           enddo
!
!
!
 6004     continue
          ee=e1st
          encode(80, 935, fmt) matter, e1st
!
          call histgc(bremc, .true.,
     *    ' *** use 12-1-9 function !', .true.,  ' ', .false.)
          call histgd(fmt, vo(1,1), -120,  0., 0.02, 1,
     *     '  v      !',15.0, 10., 'dx  ', nsampl, -2)
          ee=ee*estep
          infm=-3
           do   i=2,novlp
             encode(80,1315,fmt) ee
             if(i .eq. novlp) infm=-2
             call histg2(fmt, vo(1,i),-120, 0.,0.02,nsampl,
     *       pnt(i),infm)
             ee=ee*estep
           enddo
      endif
!
!
!
      if(pair) then
           do   i=1,novlp
              call histm(iva3(1,i),150,tcp(i)/nsampl/2)
           enddo
          ee=e1st
          encode(80,1935,fmt) matter, e1st
 1935     format(' matter=',a4,' fpair(v) at e=',g9.3,'gev!')
!
          call histgc(pairc, .true., ' ', .false., ' ', .false.)
          call histgd(fmt,iva3 ,-150, .5, 0.01,1,
     *     ' v      !', 15., 10., 'dx  ', nsampl, -2)
!
 1020     continue
          ee=ee*estep
          infm=-3
           do   i=2,novlp
              encode(80,1315,fmt) ee
              if(i .eq. novlp) infm=-2
              call histg2(fmt,iva3(1,i),-150,.5,0.01,nsampl,pnt(i),infm
              ee=ee*estep
           enddo
      endif
!
      call mmprld('enter t, if sample compton!', .false., retry, &5000)
      if(.not. retry) goto 5000
!         compton
 4900 continue
      eg1st=2.e-3
      call mmprrd('enter 1st energy(=2e-3)!',
     1               2.e-3, 0.01e-3, 1., eg1st, &4900)
      call setsv(iva4,1, 70*4,0)
      eg=eg1st
       do   ie1=1,4
      call compt(eg,t)
      tcc(ie1)=tprob
!
       do   j=1,nsampl
      call compt(eg,t)
      call compe(eg,eg1,ee1)
      call compa(cosg, cose)
      v=eg1/eg
      call hist(v, .0, 0.02, iva4(1,ie1), 70)
       enddo
!
      encode(80,5052,fmt) eg
 5052 format(' fcompt(v) at e=',g9.3,' gev!')
      call histw2(fmt,iva4(1,ie1), 70, 0., 0.02, 1, nsampl)
      eg=eg*3.1627
       enddo
!
      write(compc,'('' total x-section of compt'', 4g9.3,''!'')')
     *   (tcc(i),i=1,4)
!
!
       do   i=1,4
        call histm(iva4(1,i),  70, tcc(i)/nsampl)
       enddo
      eg=eg1st
      encode(80,5130, fmt) eg
 5130 format(' fcompt(v) at e=',g9.3,'gev!')
      call histgc(compc, .true., ' ', .false., ' ', .false.)
      call histgd(fmt,iva4(1,1), -70, 0.,  0.02, 1,
     *  ' v=eg1/eg!', 0, 0, 'no', 10000, -2)
      eg=eg*3.1627
      infm=-3
       do   i=2,4
      encode(80,1315,fmt) eg
      if(i .eq.  4   ) infm=-2
      call histg2(fmt,iva4(1,i),-70,0., 0.02,10000, pnt(i),infm)
      eg=eg*3.1627
       enddo
      call mmprld('enter t, if retry!',.false., retry, &4900)
      if(retry) goto 4900
!
!
 5000 continue
!
      if( gdmode) goto 5
      stop
      end
