!     ****************************************************************
!     *                                                              *
!     * bpgraph2:  display graph of probability function for brem and*
!     *            pair with landau effect at high energies          *
!     *                                                              *
!     ****************************************************************
!
!      execute cntl(oneday) to prepare gd load module
!              modify control card testbps1 and use it
!
!       *** note:  kkll.load needed
!
!      prepare:   fbrem in fbrem2 module ( -inc works for this end)
!
!
!        for brems 2 gev to 200 tev and for pair cre 2 to 200 tev.
!
!
!
      data e1stb/2./, e1stp/2000./
      logical retry
      dimension work(1000), x(500,13),y(500,13),iva(13),fmt(50)
      equivalence (de1, estep)
      integer linea(3)/'sld','sld', 'sld'/
      logical stop/.false./
      character*4 matter
!
      common /landuc/ x0cm, x0g, s1, alogs1, sconst
      common /landu1/e
!
!
!
!      material   z     a     rho     t0 (g/cm2)  t0 (cm)    ec
!         pb      82   207.2  11.35    5.82        0.513     6.7
!         cu      29    63.5   8.94   12.7        1.42     16+x
!         fe      26    55.85  7.86   13.7         1.74     18+x
!          w      74   183.92 19.3     6.275       .325      8.08
!
!
      call ggdin(0)
      call txtgin
      call txtgcd(1,1,1,5,84)
      call rdgin(1)
!
!        points are dense
      call cvgdi(.true.)
!
!
!
    1 continue
      de1=sqrt(10.)
      novlp=11
      novlpp=6
      e1stb=de1
      e1stp=e1stb*1000.
      call rdg('matter(= );e1stb(=3.162);e1stp(=3162);estep(=3.162);
     *novlp(=11); novlpp(=6);stop(l)!', 0, &15, &20)
      goto 1
   15 continue
      read(31,*,end=20) matter, e1stb,e1stp,estep, novlp, novlpp, stop
   20 continue
      call rdgcer(&1)
      if(stop) goto 9000
!
!                     ************
      if(matter .eq. 'pb') then
           z=82.
           a=207.2
           rho=11.35
           rx0=.909
      elseif(matter .eq. 'cu') then
           z=29.
           a=63.5
           rho=8.94
           rx0=.98
      elseif(matter .eq. 'fe') then
           z=26.
           a=55.85
           rho=7.86
           rx0=.986
      elseif(matter .eq. 'w ') then
           z=74.
           a=183.92
           rho=19.3
           rx0=.923
      else
          call txtgd(1, ' enter matter, z, a, density, corr.f.!',
     *    0)
   33     continue
          call rdg('matter(= ); z, a, density,rx0!', 0, &35, &38)
          goto 33
   35     continue
          read(31,*,end=38) matter, z, a, rho,rx0
   38     continue
          call rdgcer(&33)
          if(matter .eq. ' ') goto 33
      endif
      encode(100, 39, fmt)matter, z, a, rho, rx0
   39 format('matter=',a4,' z=',f7.2, ' a=',f7.2, ' rho=',g9.3,
     * ' rx0=',f6.3, '!')

      call txtgd(1, fmt, 0)
      call mmpatn(0,'wait-16')
      novlp=min0(novlp, 13)
      novlpp=min0(novlpp,13)
!
      call zpart(z, a, rho)
!
!
   27 continue
!
!        electron minimum energy (gev)
      e=e1stb
       do   ie=1, novlp
      vmin=1.e-5
      vmax=1. - 1.e-4
      iv=0
      vstep=10.**0.1
      v=vmin/vstep
   30 continue
      if(v .lt. .1) v=v*vstep
      if(v .ge. .1) v=v+.01
      v=amin1(vmax,v)
      iv=iv+1
      y(iv,ie)=fbrem(v) *rx0 *  v
      x(iv,ie)=v
      if(y(iv,ie) .le. 0.) iv=iv-1
      if(v .ne. vmax) goto 30
      iva(ie)=iv
      e=e*de1
       enddo
!
  105 continue
!         display graph
  110 continue
      e=e1stb
       do   ie1=1,novlp
      e1=e
      linec=1
!
       ie=ie1
      if(ie .ne. 1  ) goto 135
      encode(150,150,fmt)matter, z,a, rho, rx0, e1
  150 format(' matter=',a4, ' z=',f5.2,' a=',f6.2,' rho=',g9.3,
     *   ' rx0=',f5.3,     ' fbrem(v)*v at e=',g10.4,' gev!')
      call cvgd(fmt, x(1,ie), y(1,ie), 1,1, iva(ie),  1, iva(ie),
     *   6,6,'log', 'log', ' v=eg/ee!', ' fbrem(v)*v!', work,
     * 1000, -2)
       goto 155
!
  135 continue
      encode(80,156,fmt) e1
  156 format(' e=',g10.4,'!')
       linec=mod(linec,3)+1
      line=linea(linec)
      infm=-2
      call setbit(infm,1)
      if(ie   .eq. novlp) call rstbit(infm,1)
           call cvgd2(fmt, x(1,ie), y(1,ie),1,1,
     * iva(ie),1, iva(ie), '.', line, infm)
!
  155 continue
      e=e*de1
       enddo
      call mmprld('enter t,if retry!',.false.,retry,&105)
      if(retry) goto 110
!
!       ordinary brem
 5105 continue
!         display graph
 5110 continue
      e=e1stb
       do   ie1=1,novlp
      e1=e
      linec=1
!
       ie=ie1
      if(ie .ne. 1  ) goto 5135
      encode(150,5150,fmt) matter, z, a, rho,rx0, e1
 5150 format(' matter=',a4, ' z=',f5.2, ' a=',f6.2, ' rho=',g9.3,
     *  ' rx0=',f5.3, ' use .62 & .462 for rossi mag.',
     *                     ' fbrem(v)*v at e=',g10.4,' gev!')
      call cvgd(fmt, x(1,ie), y(1,ie), 1,1, iva(ie), 1, iva(ie),
     * 15.,  10.,'no', 'no',   ' v=eg/ee!', ' fbrem(v)*v!',
     * work, 1000, -2)
       goto 5155
!
 5135 continue
      encode(80,5156,fmt) e1
 5156 format(' e=',g10.4,'!')
       linec=mod(linec,3)+1
      line=linea(linec)
      infm=-2
      call setbit(infm,1)
      if(ie   .eq.  novlp) call rstbit(infm,1)
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
       do   ie=1,novlpp
      vmax=1. - 1.e-4
      vmin=1.-vmax
      iv=0
      v=vmin-.005
 1030 continue
      v=v+.005
      v=amin1(vmax,v)
      iv=iv+1
      y(iv,ie)=fpair(v)*rx0
      x(iv,ie)=v
      if(y(iv,ie) .le. 0.) iv=iv-1
      if(v .ne. vmax) goto 1030
      iva(ie)=iv
      e=e*de1
       enddo
!
!         display graph
 1110 continue
      e=e1stp
       do   ie =1,novlpp
      linec=1
!
      if(ie  .ne. 1  ) goto 1180
      encode(150,1150,fmt) matter, z, a, rho, rx0, e
 1150 format(' matter=',a4, ' z=',f5.2, ' a=',f6.2, ' rho=',g9.3,
     * ' rx0=',f5.3,  ' use .68 & 1 rossi mag.',
     *       ' fpair(v) at e=',g9.3,' gev!')
       call cvgd(fmt, x(1,ie), y(1,ie), 1,1, iva(ie), 1, iva(ie),
     * 15.,10.,'no', 'no',  ' v=ee/eg!', ' fpair(v)!', work,
     * 1000, -2)
      goto 1155
 1180 continue
      encode(80,156,fmt) e
      infm=-2
       linec=mod(linec,3)+1
      line=linea(linec)
      call setbit(infm,1)
      if(ie   .eq. novlpp ) call rstbit(infm,1)
       call cvgd2(fmt, x(1,ie), y(1,ie),1,1,
     * iva(ie), 1, iva(ie), '.', line, infm)
!
 1155 continue
      e=e*de1
       enddo
      call mmprld('enter t,if retry!',.false.,retry,&1110)
      if(retry) goto 1110
!
!
      goto 1
 9000 continue
      call ggdtm
      stop
      end
!
!       -inc fbrem2
       -inc fbrem2
 infm)
!
 1155 continue
      e=e*de1
       enddo
      call mmprld('enter t,if retry!',.false.,retry,&1110)
      if(retry) goto 1110
!
!
      goto 1
 9000 continue
      call ggdtm
      stop
      end
!
!       -inc fbrem2
       -inc fbrem2
