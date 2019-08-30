!     ****************************************************************
!     *                                                              *
!     * bpgraph1:  makes graph of brems and pair probability function*
!     *            at low energy where no landau effect exists       *
!     *            all other correction included                     *
!     *                                                              *
!     *            ee: upto    20gev                                 *
!     *            eg: upto    10gev                                 *
!     *                                                              *
!     ****************************************************************
!
!       execute cntl(oneday) to prepare gd load module
!       prepare:  fbrem ( -inc works for this end; fbrem1 module
!                              used)
!                 give z in dummy=zpart(z)
!                 give param by ft05f001: e1stb, e1stp, and stop
!
!
      common / bpcom /  al183z,e,ccz,emass,bcoef,fz ,z333
      data e1stb/3.e-3/, e1stp/4.e-3/
      logical retry
      dimension work(1000), x(500,13),y(500,13),iva(13),fmt(50)
      dimension vmaxa(13)
      integer linea(3)/'sld','sld', 'sld'/
      logical stop/.false./
!
      call ggdin(0)
      call txtgin
      call txtgcd(1,1,1,5,84)
      call rdgin(1)
!        points are dense
      call cvgdi(.true.)
      z=0.
!
!
    1 continue
      novlp=11
      de1=sqrt(10.)
      e1stb=1.e-3
      e1stp=de1*1.e-3
      call txtgd(1, 'enter 0 for z to stop!', 0)
      call rdg('z; e1stb(=1e-3); e1stp(=3.16e-3), estep(3.162),novlp(=
     *11)!',   0, &15, &20)
      goto 1
   15 continue
      read(31,*,end=20) z, e1stb,e1stp,de1, novlp
   20 continue
      if(z .le. 0.)   goto 9000
!
      dummy=zpart( z )
      novlp=min0(13, novlp)
!
!
!        electron minimum energy (gev)
!
      e=e1stb
       do   ie=1, novlp
        vmin=1.e-5
        vmax=vmaxv(e) - 1.e-4
        iv=0
        vstep=10.**0.1
        v=vmin/vstep
   30   continue
            if(v .le. .1 ) v=v*vstep
            if(v  .gt. .1 .and. v .le. .8)  v=v+.02
            if(v .gt. .8) v=v+.005
            v=amin1(vmax,v)
            iv=iv+1
            y(iv,ie)=fbrem(v) * v
            x(iv,ie)=v
            if(y(iv,ie) .le. 0.) iv=iv-1
        if(v .ne. vmax) goto 30
        iva(ie)=iv
        e=e*de1
       enddo
!
  105 continue
  110 continue
!         display graph
      e=e1stb
       do   ie1=1,novlp
      e1=e
      linec=1
!
      ie=ie1
      if(ie .ne. 1  ) goto 135
      encode(80,150,fmt) z, e1
  150 format(' z=',f7.2,
     * ' fbrem(v)*v at e=',g10.4, ' gev!')
      call cvgd(fmt, x(1,ie), y(1,ie), 1,1,  iva(ie), 1, iva(ie),
     *    6,6,'log', 'log', ' v=eg/ee!', ' fbrem(v)*v!', work,
     * 1000, -2)
       goto 155
!
  135 continue
      encode(80,156,fmt) e1
  156 format(' e=',g10.4,' gev!')
       linec=mod(linec,3)+1
      line=linea(linec)
      infm=-2
      call setbit(infm,1)
      if(ie   .eq. novlp ) call rstbit(infm,1)
           call cvgd2(fmt, x(1,ie), y(1,ie),1,1,
     * iva(ie), 1, iva(ie), '.', line, infm)
!
  155 continue
      e=e*de1
       enddo
      call mmprld('enter t,if retry!',.false.,retry,&105)
      if(retry) goto 110
!
!         display graph  ord. sacle brem
 5105 continue
 5110 continue
      e=e1stb
       do   ie1=1,novlp
      e1=e
      linec=1
!
      ie=ie1
      if(ie .ne. 1  ) goto 5135
      encode(150,5150,fmt) z, e1
 5150 format(' z=',f7.2, ' use .62 & .462 for x & y mag. for rossi',
     * ' fbrem(v)*v at e=',g10.4, ' gev!')
      call cvgd(fmt, x(1,ie), y(1,ie), 1,1, iva(ie), 1, iva(ie),
     *  15., 10.,'no', 'no', ' v=eg/ee!', ' fbrem(v)*v!',
     * work, 1000, -2)
       goto 5155
!
 5135 continue
      encode(80,5156,fmt) e1
 5156 format(' e=',g10.4,' gev!')
       linec=mod(linec,3)+1
      line=linea(linec)
      infm=-2
      call setbit(infm,1)
      if(ie   .eq. novlp ) call rstbit(infm,1)
           call cvgd2(fmt, x(1,ie), y(1,ie),1,1,
     * iva(ie),1, iva(ie), '.', line, infm)
!
 5155 continue
      e=e*de1
       enddo
      call mmprld('enter t,if retry!',.false.,retry,&5105)
      if(retry) goto 5110
!
!         pair
!
      e=e1stp
       do   ie=1, novlp
      vmax=vmaxv(e) - 1.e-4
      vmaxa(ie)=vmax
      vmin=1.-vmax
      iv=0
      v=vmin-.005
 1030 continue
      v=v+.005
      v=amin1(vmax,v)
      iv=iv+1
      y(iv,ie)=fpair(v)
      x(iv,ie)=v
      if(y(iv,ie) .le. 0.) iv=iv-1
      if(v .ne. vmax) goto 1030
      iva(ie)=iv
      e=e*de1
       enddo
!
 1105 continue
!         display graph
 1110 continue
      e=e1stp
       do   ie1=1,novlp
      e1=e
      linec=1
!
      ie=ie1
      if(ie .ne. 1  ) goto 1180
      encode(180,1150,fmt) z, e1, vmaxa(ie)
 1150 format(' z=',f7.2, ' use .68 & 1.05 for x & y mag. for ross.',
     *' fpair(v) at e=',g9.3,' gev. vmax=',f6.3,'!')
       call cvgd(fmt, x(1,ie), y(1,ie), 1,1, iva(ie), 1, iva(ie),
     *  15.,10.,'no', 'no', ' v=ee/eg!', ' fpair(v)!', work,
     * 1000, -2)
      goto 1155
 1180 continue
      encode(80,156,fmt) e1
      infm=-2
       linec=mod(linec,3)+1
      line=linea(linec)
      call setbit(infm,1)
      if(ie   .eq. novlp ) call rstbit(infm,1)
       call cvgd2(fmt, x(1,ie), y(1,ie),1,1,
     * iva(ie),1, iva(ie), '.', line, infm)
!
 1155 continue
      e=e*de1
       enddo
      call mmprld('enter t,if retry!',.false.,retry,&1105)
      if(retry) goto 1110
!
!
      goto 1
 9000 continue
      call ggdtm
      stop
      end
!
!       -inc fbrem1
       -inc fbrem1
      call mmprld('enter t,if retry!',.false.,retry,&1105)
      if(retry) goto 1110
!
!
      goto 1
 9000 continue
      call ggdtm
      stop
      end
!
!       -inc fbrem1
       -inc fbrem1
