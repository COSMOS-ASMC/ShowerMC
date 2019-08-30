c***********************************************************************
c***********************************************************************
c                                                                      *
c          PART 7: Event study and data analysis routines              * 
c                                                                      *
c             List of subprograms with main purpose                    * 
c        (S = subroutine, F = function, E = entry, B = block data)     *
c                                                                      *
c s jamana   to administer event study and analysis of JAM             *
c s jamanat  to administer time dependent event study and analysis     *
c s jamout1  to output phase space data (mdraw format)                 *
c s jamout2  to output phase space data (angle2 format)                *
c s jamfile  to organize opening and closing of output files           *
c s jamanafw to give time development of directed transverse flow      *
c s jammrun  to control multi run                                      *
c s jamcntpa to count time evelution of particles                      *
c s jamanl3  to analyze observables                                    *
c s jamanacl to analysize collision spectra                            *
c s jamadns1 to calculate time evolution of density with box shell     *
c s jamadns2 to calculate time evolution of density with Gaussian      *
c s jamlemt  to Lorentz transform energy-momentum tensor               *
c s jameksve to save kinetic energy to calculate temperature for box   *
c e jamekclr to clear memory of kinetic energy in jameksve             *
c s jamctmp  to calculate temperature in box calculation               *
c s jamtout  to output temperature in box calculation                  *
c s jamgout  to output ground state nucleus sampled                    *
c s jamana2  to count collision history and particle yield             *
c e jamclhst to count collisions at each collision                     *
c                                                                      *
c                                                                      *
c***********************************************************************
c***********************************************************************
c***********************************************************************

      subroutine jamana(msel)

c...Purpose: to administrate event stduy and analysis of JAM.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c...msel=0: reset statisitcs for the first run.
c...msel=1: reset statisitcs at each run.
c...msel=2: accumulate statstics for each time step during an event.
c...msel=3: accumulate  general statstics for each event as a whole.
c...msel=4: write out statistics at end of the run to files.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c...Reset statisitcs for the over all run.
c   ---------------------------------------
      if(msel.eq.0) then

c.......Collision history
        if(mstc(162).ne.0.or.mstc(165).ne.0) call jamana2(0)

c........Time evolution of directed transverse flow
        if(mstc(163).ge.1) call jamanafw(0)

c.......Time evolution of density
        if(mstc(166).ge.1) call jamadns1(0)
        if(mstc(167).ge.1) call jamadns2(0)

c........Observables.
        do i=2,5
          if(mstc(150+i).ge.1) call jamanl3(10*i)
        end do
c///////////
c        write(0,*) ' 888 mstc(156)=',mstc(156)
c//////////
        if(mstc(156).ge.1) call jamanacl(0)

c.......Temperature maru
        if(mstc(168).ge.1) call jamekclr(0)

c......Booking ground state.
        if(mstc(157).ge.1)  call jamgout(0,0.0d0,0)
 
c...Reset statisitcs at each run.
      else if(msel.eq.1) then

        if(mstc(162).eq.1.or.mstc(165).ne.0) call jamana2(1)
        if(mstc(163).ge.1) call jamanafw(1) ! d.t.flow
        if(mstc(166).ge.1) call jamadns1(1)
        if(mstc(167).ge.1) call jamadns2(1)
        if(mstc(168).ge.1) call jameksve(1)


c...Accumulate statstics for each time step during an event.
      else if(msel.eq.2) then

            if(mstc(162).eq.1.or.mstc(165).ne.0) call jamana2(2)
            if(mstc(163).ge.1) call jamanafw(2) ! d.t.flow

       if(mstc(164).eq.1) call jamout1(pard(1))
       if(mstc(164).eq.2) call jamout2(pard(1))
       if(mstc(166).ge.1) call jamadns1(2)
       if(mstc(167).ge.1) call jamadns2(2)
       if(mstc(168).ge.1) call jameksve(2)

c....Statistics for every collisions and decays.
      else if(msel.eq.4) then

c       if(mstc(163).ge.1) call jamanafw(4) ! d.t.flow

c...Accumulate  general statstics for each event as a whole.
      else if(msel.eq.3) then

        do i=2,5
          if(mstc(150+i).ge.1) call jamanl3(10*i+1)
          if(mstc(150+i).ge.2) call jamanl3(10*i+2)
        end do

c...Write out statistics at end of the run to files.
      else if(msel.eq.10) then

        if(mstc(162).eq.1.or.mstc(165).ne.0) call jamana2(10)
        if(mstc(163).eq.1) call jamanafw(3) ! d.t.flow

c...Output some observables( pion multiplicity, etc...)
        do i=2,5
          if(mstc(150+i).eq.1) call jamanl3(10*i+2)
        end do

        if(mstc(157).ge.1)  call jamgout(3,0.0d0,0)

        if(mstc(156).ge.1) call jamanacl(3)
        if(mstc(166).ge.1) call jamadns1(3)
        if(mstc(167).ge.1) call jamadns2(3)

c.......Temperature maru
        if(mstc(168).ge.1) call jamtout

      else
        write(check(1),'(''msel='',i9)')msel
        call jamerrm(30,1,'(jamana:) invalid msel')
      endif

      end

c***********************************************************************

      subroutine jamanat(msel)

c...Purpose: to administer time dependent event study and analysis.
      include 'jam1.inc'
      include 'jam2.inc'
      parameter(mxs=200)
      dimension nprint(mxs),mprint(mxs)
      save nprint,mprint

      if(msel.eq.0) then

        if(mstc(5).gt.mxs) then
          call jamerrm(30,0,'(jamanat:)test particle too large mxs')
        endif
        do i=1,mstc(5)
        nprint(i)=0
        end do

      else if(msel.eq.1) then

      mpri=int(pare(1)/parc(7))
      do i=1,mstc(5)
        mprint(i)=mpri
      end do

c....Check after every collision/decay.
      else if(msel.eq.2.or.msel.eq.3) then

      isimul=mstd(8)
      iprint=int(pard(1)/parc(7))
      if(iprint.ge.nprint(isimul)) then
        pard1=pard(1)
        do i=nprint(isimul),iprint
          pard(1)=i*parc(7)
          call jamana(2)
          if(mstd(8).eq.mstc(5)) then
            if(mstc(8).ge.1) call jamdisp(6,i+1)
          endif
        end do
        nprint(isimul)=iprint+1
        pard(1)=pard1
      else
        call jamana(4)
      endif

      if(msel.eq.2) return

c...Final analsis.
      if(nprint(isimul).lt.mprint(isimul)) then
        pard1=pard(1)
        do i=nprint(isimul),mprint(isimul)
          pard(1)=i*parc(7)
          call jamana(2)
          call jamana(4)
          if(mstd(8).eq.mstc(5).and.mstc(8).ge.1) call jamdisp(6,2)
        end do
        nprint(isimul)=iprint
        pard(1)=pard1
      endif

      endif

      end

c***********************************************************************

      subroutine jamout1(atime)

c...Purpose: to out put phase space data in marudraw format.
      include 'jam1.inc'
      include 'jam2.inc'
      character chfile*80,char*8,char1*12,cfile*82
      dimension rr(3,mxv),kd(mxv)
      save id
      data id/0/

      id=id+1
      if(id.gt.1000) return
c...Write file name.
      if(id.lt.10) write(chfile,400) id
      if(id.ge.10.and.id.lt.100) write(chfile,410) id
      if(id.ge.100.and.id.lt.1000) write(chfile,420) id
C...Format statements for output file.
  400 format('DISP/JAMD.00',i1)
  410 format('DISP/JAMD.0',i2)
  420 format('DISP/JAMD.',i3)
      iunit=55
      leng=index(fname(8),' ')-1
      cfile=chfile
      if(leng.gt.1) cfile=fname(8)(1:leng)//chfile
      open(unit=iunit,file=cfile,status='unknown')


      write(iunit,'(''# time '',f10.4,'' fm/c'')')atime
c      write(iunit,800)
c 800  format('X! -25 25'/'Y! -20 20')
c     write(iunit,801)
c801  format('m 4; l 0')

      mp=0
      do 100 i=1,nv
        k1=k(1,i)
        if(k1.gt.10) goto 100    ! dead particle
        if(k1.le.0) goto 100     ! within a formation time(newly)
        dt=atime-r(4,i)
c... (some const. quarks) within a formation time
        if(dt.lt.0.0d0) goto 100

        mp=mp+1
        rr(1,mp)=r(1,i)+dt*p(1,i)/p(4,i)
        rr(2,mp)=r(2,i)+dt*p(2,i)/p(4,i)
        rr(3,mp)=r(3,i)+dt*p(3,i)/p(4,i)
        kd(mp)=i
 100  continue

      do 210 i=1,mp-1
      do 200 j=i+1,mp
        if(rr(2,i).lt.rr(2,j)) goto 200
         itmp=kd(i)
         y1=rr(1,i) 
         y2=rr(2,i) 
         y3=rr(3,i) 
         rr(1,i)=rr(1,j)
         rr(2,i)=rr(2,j)
         rr(3,i)=rr(3,j)
         kd(i)=kd(j)
         rr(1,j)=y1
         rr(2,j)=y2
         rr(3,j)=y3
         kd(j)=itmp
 200  continue
 210  continue
 
        do j=1,mp
        i=kd(j)
        if(abs(k(9,i)).eq.3) then
          if(k(2,i).eq.2212.or.k(2,i).eq.2112)then
             char='nucl'
             lch=4
             char1='ol_nucl'
          else
            char='baryon'
            lch=6
            char1='ol_baryon'
          endif
        else if(abs(k(9,i)).eq.2) then
          char='diquark'
          char1='ol_diquark'
          lch=7
        else if(abs(k(9,i)).eq.1) then
          char='quark'
          char1='ol_quark'
          lch=5
        else if(k(9,i).eq.0) then
          if(k(2,i).eq.21) then
            char='gluon'
            char1='ol_gluon'
            lch=5
          else
            char='meson'
            char1='ol_meson'
            lch=5
          endif
        else
         char='other'
         char1='ol_other'
         lch=5
        endif

        if(abs(k(7,i)).eq.1) then
          char=char(1:lch)//'0;'
          char1=char1(1:lch+3)//'0;'
        else
          char=char(1:lch)//';'
          char1=char1(1:lch+3)//';'
        endif

        write(iunit,821)char,rr(3,j),rr(1,j),rr(2,j),k(2,i)
     $     ,char1,rr(3,j),rr(1,j)
c       write(iunit,820)k(i,2),(rr(j),j=1,3),(p(j,i),j=1,3),p(5,i)

      end do
 821  format(a8,3(1x,f10.4),i6,'; ',a12,2(1x,f10.4))
 822  format(a12,3(1x,f12.4),i6,1x,i3,f8.4)
 820  format(i6,3(1x,f10.3),4(1x,e15.4))

      write(iunit,'(''c 11'')')
      write(iunit,
     $ '(''t% 0.05 0.9 1.2  0.0  -1 -1 '',f10.4,''fm/c'')')atime

      close(iunit)

      end

c***********************************************************************

      subroutine jamout2(atime)

c 97/11/29 21:36
c...Purpose: to out put phase space data in ANGEL2 format.

      include 'jam1.inc'
      include 'jam2.inc'
      character chfile*30
      dimension rr(3,mxv), ior(mxv)
      save id
      data modeout/2/   ! 1:8-scene   2:animation
      data modey  /1/   ! 1:for all y 2:y=<0 only

c...Write file name
      data id/0/
      id=id+1
      iclose=0
      if (modeout.eq.1) then
        if (id.eq.1) then
          iunit=55
          open(unit=iunit,file=fname(5),status='unknown')
        elseif (id.gt.8) then
          return
        endif

      elseif (modeout.eq.2) then
        if (id.gt.1000) return
        if (id.lt.10) then 
          write(chfile,810) id
        elseif (id.lt.100) then
          write(chfile,820) id
        elseif(id.lt.1000) then
          write(chfile,830) id
        endif
810     format('DISP/JAMD.00',i1,'.ang')
820     format('DISP/JAMD.0', i2,'.ang')
830     format('DISP/JAMD.',  i3,'.ang')

        iunit=55
        open(unit=iunit,file=chfile,status='unknown')
      endif

c...Write header for ANGEL
      if (modeout.eq.1) then
        if (id.eq.1) then 
          write(iunit,*)'P:XORG(-0.3) YORG(1.5)'
        elseif (id.eq.5) then
          write(iunit,*)'Z:XORG(-3) YORG(-1.5)'
        else
           write(iunit,*)'Z:XORG(1) NOYN'
           write(iunit,*)'P:NOYT NOYN'
        endif
      elseif (modeout.eq.2) then
        write(iunit,*)'P:XORG(-0.1) YORG(-0.1)'
      endif
      if (modeout.eq.2.or.id.eq.1.or.id.eq.5) then
        write(iunit,*)'X:$z$[fm]'
        write(iunit,*)'Y:$x$[fm]'
      else
        write(iunit,*)'X:$z$[fm]'
      endif
      if (modeout.eq.1.and.id.eq.1) then
        write(iunit,*)'P:SCAL(0.45) FORM(1) NOFR NOMS'
      elseif (modeout.eq.2) then 
        write(iunit,*)'P:SCAL(1.5)  FORM(0.8) NOFR NOMS'
      endif
      write(iunit,*)'P:CLIN(K) CLSL(E)'
      write(iunit,*)'P:XMIN(-25) XMAX(25) YMIN(-20) YMAX(20)'
      write(iunit,'(''W:$t=$'',f10.3,'' fm/$c$/X(-20) Y(15)'')') atime
      write(iunit,*)'H:Y,NE17 N X'
      write(iunit,*)'1000 1000 1000'
      write(iunit,*)'H:Y,NB17 N X'
      write(iunit,*)'1000 1000 1000'
      write(iunit,*)'H:Y,NR17 N X'
      write(iunit,*)'1000 1000 1000'
      write(iunit,*)'H:Y,NXXXG17 N X'
      write(iunit,*)'1000 1000 1000'
      write(iunit,*)'H:Y,NXXXY17 N X'
      write(iunit,*)'1000 1000 1000'

c...Transport each particle drawn
      ii=0
      do 100 i=1,nv
        if (k(1,i).le.0.or.k(1,i).gt.10) goto 100 ! dead particle
        if (modey.eq.2.and.r(2,i).gt.0)  goto 100
        dt=atime-r(4,i)
        if(dt.lt.0.0) goto 100  ! within a formation time
        ii=ii+1

        do j=1,3
          rr(j,i)=r(j,i)+dt*p(j,i)/p(4,i)
        end do 

c...Ordering along y-axis, rr(2,ior(1))<rr(2,ior(2))<...
        if (ii.eq.1) then ! 1 st particle
          ior(1)=i
        elseif (rr(2,i).lt.rr(2,ior(1))) then    ! deepest
          do l=ii-1,1,-1
            ior(l+1)=ior(l)
          enddo
          ior(1)=i
        elseif (rr(2,i).ge.rr(2,ior(ii-1))) then ! shallowest
          ior(ii)=i
        else
          do j=1,ii-1
            if (rr(2,i).ge.rr(2,ior(j)).and.
     &          rr(2,i).lt.rr(2,ior(j+1))) then
              do l=ii-1,j+1,-1
                ior(l+1)=ior(l)
              enddo
              ior(j+1)=i
              goto 100
            endif
          enddo
          write(6,*)'(jamout2:)something is wrong'
          stop
        endif

100   continue

c...Write out cordinate in order y-component
      do i=1,ii
        if (abs(k(9,ior(i))).eq.3) then
          if (k(2,ior(i)).eq.2212.or.k(2,ior(i)).eq.2112) then !nucleon
            write(iunit,*)'H:Y,NB17 N X'
          else                          ! baryon
            write(iunit,*)'H:Y,NY17 N X'
          endif
        elseif (abs(k(9,ior(i))).eq.1) then ! quark
          write(iunit,*)  'H:Y,NR17XXX N X'
        elseif (abs(k(9,ior(i))).eq.2) then ! diquark
          write(iunit,*)  'H:Y,NR17XX N X'
        elseif (k(9,ior(i)).eq.0) then
          if (k(2,ior(i)).eq.21) then
            write(iunit,*)'H:Y,NR17XXX N X' ! gluon
          else
            write(iunit,*)'H:Y,NG17XX N X'  ! meson
          endif
        else
          write(iunit,*)'H:Y,NG17XXX N X'
        endif

        write(iunit,821)rr(1,ior(i)),rr(2,ior(i)),rr(3,ior(i))
      enddo
821   format(3(1x,f12.4))

      if ((modeout.eq.1.and.id.eq.8).or.modeout.eq.2) close(iunit)

      end

c***********************************************************************        

      subroutine jamfile(key,ifile,ctag)

c...Purpose: to organize opening and closing of output files.
c...key: =0: open, 1:close

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      character chfile*80,chfile1*81,ctag*1,cfile*82

c...Error checks.
      if(key.lt.0.or.key.gt.1) goto 110
      if(ifile.lt.100.or.ifile.gt.999) goto 130

c...Open file.
      if(key.eq.0) then

c...Write file name.
      if(ifile.eq.100) then
        ic=1
        write(chfile,100)
100     format ('JAMRUN.DAT')
      else
        if(ctag(1:1).eq.' ') then
          ic=1
          write(chfile,105) ifile
105       format ('JAM',i3,'.DAT')
        else
          ic=2
          if(ctag(1:1).eq.'a') then
            write(chfile1,106) ifile
106         format ('JAM',i3,'a.DAT')
          else if(ctag(1:1).eq.'b') then
            write(chfile1,107) ifile
107         format ('JAM',i3,'b.DAT')
          else if(ctag(1:1).eq.'c') then
            write(chfile1,108) ifile
108         format ('JAM',i3,'c.DAT')
          else if(ctag(1:1).eq.'d') then
            write(chfile1,109) ifile
109         format ('JAM',i3,'d.DAT')
          endif
        endif
      endif

        iunit=mstc(36)
        if(ic.eq.1) then
          leng=index(fname(8),' ')-1
          cfile=chfile
          if(leng.gt.1) cfile=fname(8)(1:leng)//chfile
          open(unit=iunit,file=cfile,status='unknown')
        else
          leng=index(fname(8),' ')-1
          cfile=chfile1
          if(leng.gt.1) cfile=fname(8)(1:leng)//chfile1
          open(unit=iunit,file=cfile,status='unknown')
        endif

c...Cloes file.
      else if(key.eq.1) then
        iunit=mstc(36)
        close(iunit)
      endif

      return

c...Error exit: no action taken.
110   continue  
      call jamerrm(11,0,'(jamfile): unknown key')
      return
130   continue
      call jamerrm(11,0,'(jamfile): undefined file number')

      end


c***********************************************************************

      subroutine jammrun(msel,mrun)

c...Purrpose: to control multi run.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      mrun=0
      leng=index(fname(4),' ')-1
c...Initialize control muli run, open file
      if(msel.eq.0) then

        if(mstc(39).ge.1) then
          open (mstc(39),file=fname(4)(1:leng),status='unknown')
          write(mstc(39),*)' 0 *'
          write(mstc(39),*)' 1'
          close(mstc(39))
        endif

      else if(msel.eq.1) then
        open (mstc(39),file=fname(4)(1:leng),status='old')
        read (mstc(39),*)
        read (mstc(39),*) icont
        close(mstc(39))

c....Stop execution by the use.
        if(icont.eq.0) then
           mrun=1
           call jamfin
           stop
        endif

c...Output the present event number.
        open (mstc(39),file=fname(4)(1:leng),status='old')
        write(mstc(39),*) ' ',mstd(21),' *'
        write(mstc(39),*) ' 1'
        close(mstc(39))

      endif

      end

c***********************************************************************

      subroutine jamanafw(msel)

c...Purpose: to give time development of directed transverse flow.
c...directed flow of hadrons.
c...<v1> hadrons
c...<v2> hadrons
c...<px^2> hadrons
c...<py^2> hadrons

      include 'jam1.inc'
      include 'jam2.inc'
      parameter(mtime=100,ma=4)
      dimension ntdet(0:mtime),dtm(0:mtime,ma),v1(0:mtime,ma),
     $ v2(0:mtime,ma),pxs(0:mtime,ma),pys(0:mtime,ma)
      dimension dtmx(0:mtime,ma),v1x(0:mtime,ma),
     $ v2x(0:mtime,ma),pxsx(0:mtime,ma),pysx(0:mtime,ma)

      dimension ndet(0:mtime,ma),mdet(0:mtime,ma)
      dimension dtm0(ma),v10(ma),v20(ma),pxs0(ma),pys0(ma),np(ma),mp(ma)
      double precision dbetpr
      save dtm,v1,v2,pxs,pys,ycut,yproj
      save dtmx,v1x,v2x,pxsx,pysx
      save ictime,maxtim,ndet,mdet,ntdet

c...Reset counters of time evolution of flow.
      if(msel.eq.0) then

        do i = 0, mtime
          do j=1,ma
           dtm(i,j) = 0.0d0
           v1(i,j) = 0.0d0
           v2(i,j) = 0.0d0
           pxs(i,j) = 0.0d0
           pys(i,j) = 0.0d0

           dtmx(i,j) = 0.0d0
           v1x(i,j) = 0.0d0
           v2x(i,j) = 0.0d0
           pxsx(i,j) = 0.0d0
           pysx(i,j) = 0.0d0

           ndet(i,j)=0
           mdet(i,j)=0
          end do
          ntdet(i)=0
        end do
        dbetpr=pard(35)
c       dbetta=pard(45)
        yproj=0.5d0*log((1.0d0+dbetpr)/(1.0d0-dbetpr))
c       ycut=0.3*yproj
        ycut=0.8d0*yproj

      else if(msel.eq.1) then

        ictime = 0
        maxtim=0

c...Calculate the directed transverse momenta.
      else if(msel.eq.2) then
c     else if(msel.eq.4) then

        ictime=int(pard(1)/parc(7))
        if(ictime.gt.mtime) return
        ntdet(ictime)=ntdet(ictime)+1
        maxtim=max(maxtim,ictime)

c...Determine reaction plane.
c       vecr=sqrt(pard(7)**2+pard(8)**2)
c       e1=pard(7)/vecr
c       e2=pard(8)/vecr

        do i=1,ma
        dtm0(i)=0.0d0
        v10(i)=0.0d0
        v20(i)=0.0d0
        pxs0(i)=0.0d0
        pys0(i)=0.0d0
        np(i)=0
        mp(i)=0
        end do

c....Loop over particles.
        do 300 i = 1, nv

c         if(k(1,i).le.0.or.k(1,i).gt.10) goto 300
          if(k(1,i).eq.0.or.k(1,i).gt.10) goto 300
          kf=k(2,i)
          kfl2=mod(abs(kf)/10,10)
          kfl4=mod(abs(kf)/1000,10)
          if(abs(kf).le.100.or.kfl2.eq.0) goto 300

          rap=0.5d0*log((p(4,i) + p(3,i) )/( p(4,i) - p(3,i)))
          sgn=sign(1.0d0,rap)
          px=p(1,i)
          py=p(2,i)
          pt=sqrt(px**2+py**2)
          if(pt.lt.0.03d0) goto 300
          cos1=px/pt
          cos2=(px/pt)**2-(py/pt)**2

c....Const. quarks.
          if((pard(1).lt.r(5,i)).or.(pard(1).lt.r(4,i))) then
            np(4)=np(4)+1
            dtm0(4)=dtm0(4)+sgn*px
            v10(4)=v10(4)+cos1*sgn
            if(abs(rap).le.ycut) then
              mp(4)=mp(4)+1
              v20(4)=v20(4)+cos2
              pxs0(4)=pxs0(4)+px**2
              pys0(4)=pys0(4)+py**2
            endif
            goto 300
          endif

c......All hadrons.
          np(1)=np(1)+1
c         dtm0(1)=dtm0(1)+sgn*(px*e1+py*e2)
          dtm0(1)=dtm0(1)+sgn*px
          v10(1)=v10(1)+cos1*sgn

c...Rapidity cut.
          if(abs(rap).le.ycut) then
            mp(1)=mp(1)+1
            v20(1)=v20(1)+cos2
            pxs0(1)=pxs0(1)+px**2
            pys0(1)=pys0(1)+py**2
          endif

c.....All baryons.
          if(k(9,i).eq.3) then
            j=2
          else
            j=3
          endif

          np(j)=np(j)+1
          dtm0(j)=dtm0(j)+sgn*px
          v10(j)=v10(j)+cos1*sgn
          if(abs(rap).le.ycut) then
            mp(j)=mp(j)+1
            v20(j)=v20(j)+cos2
            pxs0(j)=pxs0(j)+px**2
            pys0(j)=pys0(j)+py**2
          endif

c....End loop over particles.
 300    continue
 
        do j=1,ma
          if(np(j).ge.1) then
             ndet(ictime,j)=ndet(ictime,j)+np(j)
             dtm(ictime,j)=dtm(ictime,j)+dtm0(j)/dble(np(j))
             v1(ictime,j)=v1(ictime,j)+v10(j)/dble(np(j))
             dtmx(ictime,j)=dtmx(ictime,j)+dtm0(j)
             v1x(ictime,j)=v1x(ictime,j)+v10(j)
          endif
          if(mp(j).ge.1) then
            mdet(ictime,j)=mdet(ictime,j)+mp(j)
            v2(ictime,j)=v2(ictime,j)+v20(j)/dble(mp(j))
            pxs(ictime,j)=pxs(ictime,j)+pxs0(j)/dble(mp(j))
            pys(ictime,j)=pys(ictime,j)+pys0(j)/dble(mp(j))
            v2x(ictime,j)=v2x(ictime,j)+v20(j)
            pxsx(ictime,j)=pxsx(ictime,j)+pxs0(j)
            pysx(ictime,j)=pysx(ictime,j)+pys0(j)
          endif
        end do

c...Write the directed transverse flow momenta.
      else if(msel.eq.3) then 

          wei=1.0d0/dble(mstd(21)*mstc(5))
          ida=mstc(36)

c....Print out hadrons.
          do j=1,ma

          if(j.eq.1) then
             call jamfile(0,163,'a')
             write(ida,800)yproj,ycut
          else if(j.eq.2) then
             call jamfile(0,163,'b')
             write(ida,810)yproj,ycut
          else if(j.eq.3) then
             call jamfile(0,163,'c')
             write(ida,820)yproj,ycut
          else if(j.eq.4) then
             call jamfile(0,163,'d')
             write(ida,830)yproj,ycut
          endif

          do i=0,maxtim
            wei1=0.0d0
            wei2=0.0d0
            if(ndet(i,j).ge.1) wei1=1.0d0/ndet(i,j)
            if(mdet(i,j).ge.1) wei2=1.0d0/mdet(i,j)

            aaa=pxs(i,j)+pys(i,j)
            alp=0.0d0
            if(aaa.gt.1d-6) then
              alp=(pxs(i,j)-pys(i,j))/aaa
            endif

            aaa=pxsx(i,j)+pysx(i,j)
            alp1=0.0d0
            if(aaa.gt.1d-6) then
              alp1=(pxsx(i,j)-pysx(i,j))/aaa
            endif

            write(ida,901)i*parc(7),wei*dtm(i,j),wei*v1(i,j)
     $ ,wei*v2(i,j),wei*pxs(i,j),wei*pys(i,j)
     & ,alp,ndet(i,j),mdet(i,j),ntdet(i)
     $ ,dtmx(i,j)*wei1,v1x(i,j)*wei1,v2x(i,j)*wei2,alp1
          end do

          if(j.eq.1) call jamfile(1,163,'a')
          if(j.eq.2) call jamfile(1,163,'b')
          if(j.eq.3) call jamfile(1,163,'c')
          if(j.eq.4) call jamfile(1,163,'d')
        end do

      endif

800   format('# hadrons: yproj=',f7.3,' ycut=',f7.3,/,
     $ '# time(fm/c) p_{x}^{dir} <v1> <v2> <px^2> <py^2> alp')
810   format('# baryons: yproj=',f7.3,' ycut=',f7.3,/,
     $ '# time(fm/c) p_{x}^{dir} <v1> <v2> <px^2> <py^2> alp')
820   format('# mesons: yproj=',f7.3,' ycut=',f7.3,/,
     $ '# time(fm/c) p_{x}^{dir} <v1> <v2> <px^2> <py^2> alp')
830   format('# const.quarks: yproj=',f7.3,' ycut=',f7.3,/,
     $ '# time(fm/c) p_{x}^{dir} <v1> <v2> <px^2> <py^2> alp')
901   format(f7.2,6(1x,g12.4),3(1x,i6),4(1x,g12.4))
 
      end

c***********************************************************************
c
c   Analysis of collision history
c
c***********************************************************************

      subroutine jamana2(msel)

c...Purpose: to initiallize counters for collision history.
c...time dep. of Collisions
      include 'jam1.inc'
      include 'jam2.inc'

      parameter(mxtim=300)
      parameter(mxpart=10,mxcl=17)
      common/anldat1/nopart(mxpart,mxtim),iquark(-6:6,mxtim),
     $ iopart(10,mxtim),ncoll(0:mxcl,mxtim),lcoll(0:mxcl)
      save /anldat1/
      save ictime

c....Initialize counters.
      if(msel.eq.0) then

      do it=1,mxtim
        do i=1,mxpart
          nopart(i,it)=0
        end do
        do i = 0,mxcl
          ncoll(i,it)=0 
        end do
        do i=-6,6
          iquark(i,it)=0
        end do
        do i=1,10
          iopart(i,it)=0
        end do
      end do

c...Initilize collision counter at each run
      else if(msel.eq.1) then

      ictime = 0
      do i=0,mxcl
        lcoll(i)=0 
      end do
 
c...Count particle number and collisions.
      else if(msel.eq.2) then

        atime=pard(1)
        ictime=ictime+1
        if(ictime.gt.mxtim) return

c....Count particles.
      if(mstc(165).ne.0) then
        call jamcntpa(atime,ictime,msel)
      endif

c...Count collisions.
      if(mstc(162).ne.0) then
        do i=0,mxcl
          ncoll(i,ictime)=ncoll(i,ictime)+lcoll(i)
          lcoll(i)=0
        end do
      endif


c...Print evolution of particle and collision histroy.
      else if(msel.eq.10) then

      idp=mstc(36)
      wei=1.d0/dble(mstd(21)*mstc(5))
      ict=min(mxtim,ictime)

c....Time evolution of particles.
      if(mstc(165).ne.0) then

        call jamfile(0,165,' ')
c       write(idp,800)(nopart(j,0)*wei,j=1,mxpart)
c800    format('#time(fm/c)  N',f10.3,' Delta',f10.3,' B*',f10.3,
c    $      ' pion',f10.3,' M*',f10.3,' quark',f10.3,' gluon',
c    $        f10.3,' parton',f10.3)

        write(idp,801)
 801    format('# time(fm/c) nucl.  delta  B*   pi   M'
     $       ,' quark  gluon  parton  const.quarks other')
        do it=1,ict
          write(idp,810)it*parc(7),(nopart(j,it)*wei,j=1,mxpart)
        end do
        call jamfile(1,165,' ')
810     format(f8.4,20(1x,f10.3))

c...Print out quark content.
        call jamfile(0,165,'a')
        write(idp,802)
 802    format('# time(fm/c) d u s ratio')
        do it=1,ict
          rsd=0.0d0
          totud=iquark(-1,it)+iquark(1,it)+iquark(-2,it)+iquark(2,it)
          if(totud.gt.1d-5) then
           rsd=200*(iquark(-3,it)+iquark(3,it))/totud
          endif
          write(idp,820)it*parc(7),
     $    ((iquark(-j,it)+iquark(j,it))*wei,j=1,3),rsd
        end do
820     format(f8.4,20(1x,f10.5))
        call jamfile(1,165,'a')

c...Print out hadron content.
        call jamfile(0,165,'b')
        write(idp,803)
 803    format('# time(fm/c) pi  str   Y   xi   omega')
        do it=1,ict
          write(idp,820)it*parc(7),
     $    (iopart(j,it)*wei,j=1,6)
        end do
        call jamfile(1,165,'b')



      endif

c...Print collision histroy.
      if(mstc(162).ne.0) then
        call jamfile(0,162,' ')
        write(idp,900)
900     format('#  1) time(fm/c)'/
     $         '# 2)total'/'# 3)elastic'/'# 4)NN->NR'/'# 5)NN->RR'/
     $         '# 6)NR->NN'/'# 7)RR->NN'/'# 8)RR->NR'/
     $         '# 9)RR->RR,NR->NR'/'# 10)MB->R'/'# 11)MB->MB'/
     $         '# 12)MM->M'/'# 13)MM->MM'/
     $         '# 14)jet decay'/'# 15) resonance decay'/
     $         '# 16)BB 17) MB 18)MM 19)aBB')

        do it=1,ict
          write(idp,910)it*parc(7),(ncoll(j,it)*wei,j=0,mxcl)
        end do
        call jamfile(1,162,' ')
910     format(f8.4,20(1x,e10.4))
      endif

      endif

      return

c***********************************************************************

      entry jamclhst(ic,ichanel)

c   Purpose:  count the collision number.
c...called by coll.f, decay.f
c....0)total
c....1)elastic
c.. .2)NN->NR
c....3)NN->RR'
c....4)NR->NN
c....5)RR->NN
c....6)RR->NR'
c....7)RR->R'R'',NR->NR'
c....8)MB->R
c....9)MB->M'B'
c...10)MM->M'
c...11)MM->M'M'
c...12) jet decay
c...13) resonance decay
c...14) BB collisions
c...15) MB collisions
c...16) MM collisions
c...17) aBB collisions

c....Decay event.
      if(ic.ge.2) then
        if(ic.eq.2) lcoll(12)=lcoll(12)+1
        if(ic.eq.3) lcoll(13)=lcoll(13)+1
        return
      endif

      icltyp=mste(2)
      lcoll(0)=lcoll(0)+1
      imes1=mod(kcp(2,1)/1000,10)
      imes2=mod(kcp(2,2)/1000,10)
      if(icltyp.eq.1)lcoll(14)=lcoll(14)+1
      if(icltyp.eq.2)lcoll(15)=lcoll(15)+1
      if(icltyp.eq.3)lcoll(16)=lcoll(16)+1
      if(icltyp.eq.4)lcoll(17)=lcoll(17)+1

c...This is elastic collision.
      if(ichanel.eq.1) then
        lcoll(1)=lcoll(1)+1
        return
c....MB->R
      else if(ichanel.eq.3) then
        if(imes1.eq.0.and.imes2.eq.0) then
          lcoll(10)=lcoll(10)+1
        else
          lcoll(8)=lcoll(8)+1
        endif
        return
      endif

c...MM->M'M'
      if(imes1.eq.0.and.imes2.eq.0) then
        lcoll(11)=lcoll(11)+1
        return
c...MB->M'B'
      else if(imes1.eq.0.or.imes2.eq.0)then
        lcoll(9)=lcoll(9)+1
        return
      endif

c....Inititial particle IDs
      id1=kchg(mste(22),5)
      id2=kchg(mste(24),5)
c     if(id1.eq.id_nucl.and.pare(25).gt.parc(25)) id1=id_nucls
c     if(id1.eq.id_pi.and.pare(25).gt.parc(26)) id1=id_rho
c     if(id2.eq.id_nucl.and.pare(50).gt.parc(25)) id2=id_nucls
c     if(id2.eq.id_pi.and.pare(50).gt.parc(26)) id2=id_rho

c...Final particle IDs
      id3=0
      id4=0
      if(mste(26).gt.0) id3=kchg(mste(26),5)
      if(mste(28).gt.0) id4=kchg(mste(28),5)
c     if(id3.eq.id_nucl.and.pare(75).gt.parc(25)) id3=id_nucls
c     if(id3.eq.id_pi.and.pare(75).gt.parc(26)) id3=id_rho
c     if(id4.eq.id_nucl.and.pare(99).gt.parc(25)) id4=id_nucls
c     if(id4.eq.id_pi.and.pare(99).gt.parc(26)) id4=id_rho

c....NN -> X
      if(id1.eq.id_nucl.and.id2.eq.id_nucl) then
        if(id3.eq.id_nucl.and.id4.ne.id_nucl) then
          lcoll(2)=lcoll(2)+1
        else if(id3.ne.id_nucl.and.id4.eq.id_nucl) then
          lcoll(2)=lcoll(2)+1
        else
          lcoll(3)=lcoll(3)+1
        endif
      else 

c......NR -> X
        if(id1.eq.id_nucl.or.id2.eq.id_nucl) then
c.......NR ->NN
          if(id3.eq.id_nucl.and.id4.eq.id_nucl) then
            lcoll(4)=lcoll(4)+1
          else
            lcoll(7)=lcoll(7)+1
          endif

c......RR -> X
        else

c.......RR ->NN
          if(id3.eq.id_nucl.and.id4.eq.id_nucl) then
            lcoll(5)=lcoll(5)+1
c.......RR ->RR
          else if(id3.ne.id_nucl.and.id4.ne.id_nucl) then
            lcoll(7)=lcoll(7)+1
c.......RR ->RN
          else
            lcoll(6)=lcoll(6)+1
           endif
        endif
  
      endif

      end

c***********************************************************************

      subroutine jamcntpa(atime,ict,msel)

c...Count time evelution of parton/particle.
      include 'jam1.inc'
      include 'jam2.inc'
      parameter(mxtim=300)
      parameter(mxpart=10,mxcl=17)
      common/anldat1/nopart(mxpart,mxtim),iquark(-6:6,mxtim),
     $ iopart(10,mxtim),ncoll(0:mxcl,mxtim),lcoll(0:mxcl)
      save /anldat1/

c..1)nucl. 2)delta(1232) 3)Bs* 4)pi 5)M 6)quark 7)gluon 8)parton

c     WRITE(30,*)'time nv',atime,nv
      nq=0

      do 100 i=1,nv
c       WRITE(30,800)k(1,i),k(2,i),r(4,i),r(5,i)
c800    format(i3,1x,i9,2(1x,g10.3))

       if(k(1,i).ge.11) goto 100
       if(r(4,i).gt.atime.and.msel.eq.2) goto 100

c...number of const. quarks.
       if(r(5,i).gt.atime.and.msel.eq.2) then
         iq=mod(abs(k(1,i))/10,10)
         if(iq.eq.3) iq=2
         nopart(9,ict)=nopart(9,ict)+iq
c          nq=nq+iq
c          if(iq.ne.1.and.iq.ne.2)write(30,*)'iq=',iq,k(1,i),k(2,i)
c          if(k(1,i).ge.0)write(30,*)'iq=',iq,k(1,i),k(2,i)
         goto 100
       endif

       kf=k(2,i)
       em=p(5,i)
       kfa=abs(kf)
       kc=jamcomp(kf)
       id=kchg(kc,5)
       k9=k(9,i)
       k9a=abs(k9)

c      if(k9a.eq.3) nq=nq+3
c      if(k9a.eq.0) nq=nq+1

c...nopart(1,):nucl.
c...nopart(2,):delta
c...nopart(3,):B*
c...nopart(4,):pion
c...nopart(5,):all meosons
c...nopart(6,):quark
c...nopart(7,):gluon
c...nopart(8,):patons
c...nopart(9,): const. quarks
c...nopart(10,):others

       if(kfa.le.100) then
         if(kfa.lt.10) then
           nopart(6,ict)=nopart(6,ict)+1
           nopart(8,ict)=nopart(8,ict)+1
         else if(kf.eq.21) then
           nopart(7,ict)=nopart(7,ict)+1
           nopart(8,ict)=nopart(8,ict)+1
         else
           nopart(10,ict)=nopart(10,ict)+1
         endif
       else if(mod(kfa/10,10).eq.0) then ! diquark
         nopart(6,ict)=nopart(6,ict)+2
         nopart(8,ict)=nopart(8,ict)+2

       else if(k9a.eq.3) then

       if(id.eq.id_nucl) then
         nopart(1,ict)=nopart(1,ict)+1
       else
         if(id.eq.id_delt) then
           nopart(2,ict)=nopart(2,ict)+1
         else if(em .gt.1.232) then
           nopart(3,ict)=nopart(3,ict)+1
         endif
       endif

       else if(k9a.eq.0) then
         if(id.eq.id_pi) then
           nopart(4,ict)=nopart(4,ict)+1
         else if(em.gt.1.232) then
           nopart(5,ict)=nopart(5,ict)+1
         endif
       else
         nopart(10,ict)=nopart(10,ict)+1
       endif

c...Rapidiy cut.
      y=0.5d0*log( max(p(4,i)+p(3,i),1.d-8)/max(p(4,i)-p(3,i),1.d-8) )
      if(abs(y).gt.1.0) goto 100

c....Get quark contenet of the hadron.
      call attflv2(kf,ifla,iflb,iflc)
      if(ifla.ne.0) iquark(ifla,ict)=iquark(ifla,ict)+1
      if(iflb.ne.0) iquark(iflb,ict)=iquark(iflb,ict)+1
      if(iflc.ne.0) iquark(iflc,ict)=iquark(iflc,ict)+1

c...Count particles.
       jj=0
       if(id.eq.id_pi)     jj=1
       if(id.eq.id_light1) jj=1
       if(id.eq.id_light0) jj=1
       if(id.eq.id_str)    jj=2
       if(id.eq.id_lamb)   jj=3
       if(id.eq.id_lambs)  jj=3
       if(id.eq.id_sigm)   jj=3
       if(id.eq.id_sigms)  jj=3
       if(id.eq.id_xi)     jj=4
       if(id.eq.id_xis)    jj=4
       if(id.eq.id_omega)  jj=5
       if(jj.ne.0) iopart(jj,ict)=iopart(jj,ict)+1 

c      if(id.eq.id_pi)     iopart(1,ict)=iopart(1,ict)+1 
c      if(id.eq.id_light1) iopart(1,ict)=iopart(1,ict)+1 
c      if(id.eq.id_light0) iopart(1,ict)=iopart(1,ict)+1 
c      if(id.eq.id_str)    iopart(2,ict)=iopart(2,ict)+1 
c      if(id.eq.id_lamb)   iopart(3,ict)=iopart(3,ict)+1 
c      if(id.eq.id_lambs)  iopart(3,ict)=iopart(3,ict)+1 
c      if(id.eq.id_sigm)   iopart(3,ict)=iopart(3,ict)+1 
c      if(id.eq.id_sigms)  iopart(3,ict)=iopart(3,ict)+1 
c      if(id.eq.id_xi)     iopart(4,ict)=iopart(4,ict)+1 
c      if(id.eq.id_xis)    iopart(4,ict)=iopart(4,ict)+1 
c      if(id.eq.id_omega)  iopart(5,ict)=iopart(5,ict)+1 
c      if(kf.eq.321)  iopart(6,ict)=iopart(6,ict)+1 


100   continue

c      write(30,*)nq/3.0d0

      end

************************************************************************
      subroutine jamanl3(msel)

c...Purpose: to analyze observables.

      include 'jam1.inc'
      include 'jam2.inc'

c....Rapidity distribution
      parameter(mpar=2)
      parameter(ymax=2.0d0,ymin=-2.0d0,my=15)
      dimension ny(mpar,my),pxave(mpar,2,my),ptave(mpar,2,my)
      save pxave,ptave,ny
      save  ylab, ycmpr,ypr,dy,yproj

c...Reset statistics on rapidity distributions dn/dy.
      if(msel.eq.20) then
      else if(msel.eq.21) then
      else if(msel.eq.22) then

c...Reset statistics on dn/dp_t distributions.
      else if(msel.eq.30) then
      else if(msel.eq.31) then
      else if(msel.eq.32) then

c...Reset statistics on pion multiplicity.
      else if(msel.eq.40) then

c...Count pion multiplicity
      else if(msel.eq.41) then

c...Print multiplicity.
      else if(msel.eq.42) then

c-----------------------------------------
c...Reset statisitcs on transeverse flow.
c-----------------------------------------
      else if(msel.eq.50) then

        dy=(ymax-ymin)/my
        ylab=pard(17)
        betacm=pard(5)
        elab=pard(14)
        plab=pard(15)
        ylab = 0.5d0 * log( ( 1.0d0 + betacm ) / ( 1.0d0 - betacm ) )
c       plab = sqrt( elab * elab + 2. * em_nuc * elab )
        betapr = plab / ( elab + pard(44) )
        ypr = 0.5d0 * log( ( 1.0d0 + betapr ) / ( 1.0d0 - betapr ) )
        ycmpr = ypr - ylab
        dbetpr=pard(35)
        yproj=0.5d0*log((1.0d0+dbetpr)/(1.0d0-dbetpr))

        do i=1,my
          do j=1,mpar
          ny(j,i)=0
          do l=1,2
          pxave(j,l,i)=0.0d0
          ptave(j,l,i)=0
          end do
          end do
        end do

c...Accumulate statisitcs on <P_x(y)> and <P_T(y)>
      else if(msel.eq.51) then

c       b=pard(2)
c       xteven = xteven + pard(2)
c...Determine reaction plane.
c       vecr=sqrt(pard(7)**2+pard(8)**2)
c       e1=pard(7)/vecr
c       e2=pard(8)/vecr

c.....Loop over particles
        do 222 i = 1, nv

          if(k(1,i).gt.10) goto 222
          kf=k(2,i)
          rap=0.5d0*log( max(p(4,i)+p(3,i),1.d-8)/
     $                      max(p(4,i)-p(3,i),1.d-8) )

          pt=max(sqrt(p(1,i)**2+p(2,i)**2),1.d-8)
          iy=int((rap/yproj-ymin)/dy)+1
c         iy=int((rap-ymin)/dy)+1
          if(iy.gt.my.or.iy.le.0) goto 222

c         px=p(1,i)*e1+p(2,i)*e2 
          px=p(1,i)

          if(kf.eq.2112.or.kf.eq.2212) then
            ny(1,iy)=ny(1,iy)+1
            pxave(1,1,iy)=pxave(1,1,iy)+px
            ptave(1,1,iy)=ptave(1,1,iy)+pt
          else if(abs(kf).eq.211.or.kf.eq.111) then
            ny(2,iy)=ny(2,iy)+1
            pxave(2,1,iy)=pxave(2,1,iy)+px
            ptave(2,1,iy)=ptave(2,1,iy)+pt
          end if

222     continue

c...Print <P_x>.
      else if(msel.eq.52) then

        idb=mstc(36)

       do j=1,mpar
         if(j.eq.1) then
           call jamfile(0,155,'a')
           write(idb,'(''# nuleon flow <p_x> <p_t>'')')
         else if(j.eq.2) then
            call jamfile(0,155,'b')
            write(idb,'(''# pion flow <p_x> <p_t>'')')
         endif
         write(idb,'(''# ypr ycmpr ylab dy'',5(1x,f7.3))')
     $               ypr,ycmpr,ylab,yproj,dy

        do i = 1, my

          yrap=(i-0.5d0)*dy+ymin

          if(ny(j,i).ge.2) then
            pxt=pxave(j,1,i)/ptave(j,1,i)
            p_x=pxave(j,1,i)/ny(j,i)
            p_t=ptave(j,1,i)/ny(j,i)
            aaa=ny(j,i)*(ny(j,i)-1)
          else
            p_x=0.0d0
            p_t=0.0d0
            pxt=0.0d0
          end if
          write(idb,'(f7.3,3(1x,1pe11.3),i9)')
     $                yrap,p_x,p_t,pxt,ny(j,i)
         end do

         if(j.eq.1) call jamfile(1,155,'a')
         if(j.eq.2) call jamfile(1,155,'b')
        end do

      endif

      end

c***********************************************************************

      subroutine jamanacl(msel)

c...Analysize collision spectra.
      include 'jam1.inc'
      include 'jam2.inc'

      parameter(ns=100,mxk=15)
      dimension scoll(mxk,ns)
      dimension kfl(4,3)
      save smin,smax,ds,scoll

c...scoll(1,): BB collision
c...scoll(2,): MB collision
c...scoll(3,): MM collision
c...scoll(4,): diquark BB collision
c...scoll(5,): quark BB collision

c...Initialization.
      if(msel.eq.0) then
        if(mstc(156).eq.0) return
        smin=0.0d0
        smax=max(2.0d0*pard(16),2.0d0)
        ds=(smax-smin)/ns
        do i=1,ns
         do j=1,mxk
          scoll(j,i)=0.0d0
         end do
        end do
        return

      else if(msel.eq.2) then

      ichanel=mste(1)
      icltyp=mste(2)

      if(ichanel.eq.-1) then
        mstd(51)=mstd(51)+1
        return
      else if(ichanel.eq.-5.or.ichanel.eq.-6) then
        mstd(54)=mstd(54)+1
        return
      else if(ichanel.le.0) then
        return
c...Elastic.
      else if(ichanel.eq.1)  then
        mstd(41)=mstd(41)+1
c...Absorption
      else if(ichanel.eq.3)  then
        mstd(43)=mstd(43)+1
      endif

c...Inelastic.
      if(ichanel.ge.2) then
        mstd(42)=mstd(42)+1
c...Strangeness source.
        kf1=kcp(2,1)
        kf2=kcp(2,2)
        kf3=kcp(2,3)
        kf4=kcp(2,4)
        call attflv2(kf1,kfl(1,1),kfl(1,2),kfl(1,3))
        call attflv2(kf2,kfl(2,1),kfl(2,2),kfl(2,3))
        call attflv2(kf3,kfl(3,1),kfl(3,2),kfl(3,3))
        call attflv2(kf4,kfl(4,1),kfl(4,2),kfl(4,3))
        ist1=0
        ist2=0
        do i1=1,4
        do j1=1,3
        if(j1.le.2.and.abs(kfl(i1,j1)).eq.3) ist1=ist1+1 
        if(j1.ge.3.and.abs(kfl(i1,j1)).eq.3) ist2=ist2+1 
        end do
        end do
        if(ist2.gt.ist3) then
          mstd(61)=mstd(61)+1
          if(ks01.ge.1.and.ks02.ge.1) mstd(62)=mstd(62)+1
        endif
      endif

      ks01=kcp(1,1)
      ks02=kcp(1,2)
      if(ks01.eq.4.and.ks02.eq.4) then
        mstd(49)=mstd(49)+1 ! parton-parton collision
      else if(ks01.eq.4.or.ks02.eq.4) then
        mstd(48)=mstd(48)+1 ! parton-hadron collision
      else
        if(icltyp.eq.1) mstd(44)=mstd(44)+1 ! BB collision
        if(icltyp.eq.2) mstd(45)=mstd(45)+1 ! MB collision
        if(icltyp.eq.3) mstd(46)=mstd(46)+1 ! MM collision
        if(icltyp.eq.4) mstd(47)=mstd(47)+1 ! antiBB collision
        if(icltyp.eq.5) mstd(48)=mstd(48)+1 ! parton-hadron collision
        if(icltyp.eq.6) mstd(49)=mstd(49)+1 ! parton-parton collision
      endif

      if(mstc(156).eq.0) return

      ks11=mod(abs(ks01)/10,10)
      ks12=mod(abs(ks02)/10,10)

c....Count BB collisions.
      if(ichanel.ge.1.and.ichanel.le.5) then

        srt=pare(2)
        ix=(srt-smin)/ds
        if(ix.le.0.or.ix.gt.ns) return

c.......Count const. quark collision.
        if(ks01.ge.1.and.ks02.ge.1) then
          scoll(13,ix)=scoll(13,ix)+1d0/ds
        else
          scoll(14,ix)=scoll(14,ix)+1d0/ds
          if(ks11.eq.1.or.ks12.eq.1) scoll(15,ix)=scoll(15,ix)+1d0/ds
          if(ks11.eq.3.or.ks12.eq.3) scoll(15,ix)=scoll(15,ix)+1d0/ds
        endif

c...all BB collision.
        if(icltyp.eq.1.or.icltyp.eq.4) then
          if(icltyp.eq.1) scoll(1,ix)=scoll(1,ix)+1d0/ds

c.......Count const. quark collision for BB.
          if(ks01.ge.1.and.ks02.ge.1) then
            scoll(4,ix)=scoll(4,ix)+1d0/ds
          else
            scoll(5,ix)=scoll(5,ix)+1d0/ds
c           if(ks11.eq.2.or.ks12.eq.2) scoll(5,ix)=scoll(5,ix)+1d0/ds
            if(ks11.eq.1.or.ks12.eq.1) scoll(6,ix)=scoll(6,ix)+1d0/ds
            if(ks11.eq.3.or.ks12.eq.3) scoll(6,ix)=scoll(6,ix)+1d0/ds
          endif

c...all MB coll.
        else if(icltyp.eq.2) then
           scoll(2,ix)=scoll(2,ix)+1d0/ds

c.......Count const. quark collision for MB.
          if(ks01.ge.1.and.ks02.ge.1) then
            scoll(7,ix)=scoll(7,ix)+1d0/ds
          else
            scoll(8,ix)=scoll(8,ix)+1d0/ds
            if(ks11.eq.1.or.ks12.eq.1) scoll(9,ix)=scoll(9,ix)+1d0/ds
            if(ks11.eq.3.or.ks12.eq.3) scoll(9,ix)=scoll(9,ix)+1d0/ds
          endif

c...all MM coll.
        else if(icltyp.eq.3) then 

           scoll(3,ix)=scoll(3,ix)+1d0/ds

c.......Count const. quark collision for MB.
          if(ks01.ge.1.and.ks02.ge.1) then
            scoll(10,ix)=scoll(10,ix)+1d0/ds
          else
            scoll(11,ix)=scoll(11,ix)+1d0/ds
            if(ks11.eq.1.or.ks12.eq.1) scoll(12,ix)=scoll(12,ix)+1d0/ds
            if(ks11.eq.3.or.ks12.eq.3) scoll(12,ix)=scoll(12,ix)+1d0/ds
          endif

        endif
      endif

c...Output
      else if(msel.eq.3) then
        idb=mstc(36)
        call jamfile(0,156,' ')
        write(idb,800)
        write(idb,801)
        do ix=1,ns
          xx=smin+ds*(ix-0.5d0)
          write(idb,810)xx,(scoll(j,ix)/dble(mstd(21)),j=1,mxk)
        end do
        call jamfile(1,156,' ')

      endif

800   format('# energy distribution of collisions')
801   format('# srt(GeV) BB MB MM  BB(full B) diquark(BB)  quark(BB)')
810   format(f9.4,20(1x,f9.4))

      end

c***********************************************************************

      subroutine jamadns1(msel)

c...Purpose: calculation of density from spatial distribution of
c... testparticles and average momentum in spacial cell
c
c=======================================================================
      include 'jam1.inc'
      include 'jam2.inc'
       
      parameter(maxt=100,mab=14)
      parameter(rho0=0.168d0)
      parameter(xmax=5.0d0,ymax=5.0d0,zmax=5.0d0,vol1=xmax*ymax*zmax)
      parameter(zz=0.5d0)

      dimension a(2,mab)
      dimension current(2,4,maxt),burrent(2,4,maxt),tensor(2,4,4,maxt)
      dimension g(4,4),cu(4),cub(4),ten(4,4),cu1(4),cub1(4),ten1(4,4)
      dimension pcel(2,4,maxt),pcm(4)
      dimension tempt1(2),tempz1(2),tempt2(2),tempz2(2),temp(2)
      save kt,maxtim,current,burrent,tensor
      data g/ -1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,1/

c...istime=1: c.s.m. time
c...istime=2: proper time
      data istime/1/

c...Reset counters
      if(msel.eq.0) then
      maxtim=0
      do it = 1,maxt
        do ie=1,2
          do i=1,4
          current(ie,i,it)=0.0d0
          burrent(ie,i,it)=0.0d0
          pcel(ie,i,it)=0.0d0
          do j=1,4
           tensor(ie,i,j,it)=0.0d0
          end do
          end do
        end do
      end do

      else if(msel.eq.1) then
        kt=0

c...Calculate current and energy-momentum tensor.
      else if(msel.eq.2) then

      time=pard(1)
      if(istime.eq.1) then
        kt=kt+1
        if(kt.gt.maxt) return
        maxtim=max(maxtim,kt)
      endif


      r1=max(pard(40),1.0d0)
      r2=max(pard(50),1.0d0)
      rr2=(min(r1,r2))**2
      vol2=1.d0/(4*paru(1)*rr2*zz)

c---------------------------------------------------------------------
c....Loop over all particles
      do i=1,nv
 
        k1=k(1,i)
        if(k1.gt.10) goto 100   ! dead particle
c       if(abs(k(7,i)).eq.1) goto 100 ! not yet interaction
        if(r(5,i).gt.time) goto 100   ! not formed

        dt=time-r(4,i)
        if(dt.lt.0.0d0) goto 100

        x1=r(1,i)+dt*p(1,i)/p(4,i)
        y1=r(2,i)+dt*p(2,i)/p(4,i)
        z1=r(3,i)+dt*p(3,i)/p(4,i)
        bar=kfprop(k(2,i),2)/3.0d0

c...Proper time
        if(istime.eq.2) then
          tau2=time**2-z1**2
          if(tau2.le.0.0d0) goto 100
          tau=sqrt(tau2)
          kt=int(tau/parc(7))
          if(kt.le.0.or.kt.gt.maxt) goto 100
          maxtim=max(maxtim,kt)
        endif

        do ie=1,2
          ihit=1
          if(ie.eq.1) then
            if(abs(x1).gt.xmax/2) ihit=0
            if(abs(y1).gt.ymax/2) ihit=0
            if(abs(z1).gt.zmax/2) ihit=0
            vol=1.d0/vol1
          else
            if(abs(x1**2+y1**2).gt.rr2.or.abs(z1).gt.zz/2) ihit=0
            vol=vol2
          endif

          if(ihit.eq.1) then
          do im=1,4
            pcel(ie,im,kt)=pcel(ie,im,kt)+p(im,i)
            current(ie,im,kt)=current(ie,im,kt)+p(im,i)/p(4,i)*vol
            burrent(ie,im,kt)=burrent(ie,im,kt)+p(im,i)/p(4,i)*bar*vol
          do in=1,4
            tensor(ie,im,in,kt)=tensor(ie,im,in,kt)
     $                                +p(im,i)*p(in,i)/p(4,i)*vol
          end do
          end do
          endif
        end do

 100  continue
      end do
c---------------------------------------------------------------------

c...Output
      else if(msel.eq.3) then

        wei=1.0d0/dble(mstd(21)*mstc(5))
        ida=mstc(36)

      do ie=1,2     !( Loop over volume

        if(ie.eq.1) then
          call jamfile(0,166,'a')
          write(ida,800)
        else if(ie.eq.2) then
          call jamfile(0,166,'b')
          write(ida,801)
        endif

      do it=1,maxtim  !( Loop over time step

        do i=1,mab
         a(ie,i)=0.0d0
        end do

      do i=1,4
       cu(i)=current(ie,i,it)*wei
       cub(i)=burrent(ie,i,it)*wei
       pcm(i)=pcel(ie,i,it)
      do j=1,4
       ten(i,j)=tensor(ie,i,j,it)*wei
      end do
      end do

c...Lorentz invariant Density
      cc=cu(4)**2-(cu(1)**2+cu(2)**2+cu(3)**2)
      if(cc.gt.0.0d0) then
        a(ie,11)=sqrt(cc)
      else
c       write(6,*)'current<0',cc,(cu(j),j=1,4)
        goto 200
      endif

c...Lorentz invariant Baryon density
      bnorm=cub(4)**2 -( cub(1)**2+cub(2)**2+cub(3)**2 )
      if(bnorm.gt.0.0d0) then
        a(ie,12)=sqrt(bnorm)
      else
        a(ie,12)=0.0d0
c       write(6,*)'baryon current<0',bnorm
      endif

c...Lorentz invariant pressure and energy density
      a(ie,13)=0.0d0
      a(ie,14)=0.0d0
      do im=1,4
      do in=1,4
       tmp=g(im,in)*cu(im)*cu(in)/cc
       a(ie,13)=a(ie,13)+ten(im,in)*tmp                ! energy density
       a(ie,14)=a(ie,14)-1.d0/3.d0*ten(im,in)*(g(im,in)-tmp) ! pressure
      end do
      end do

c...Local frame 1
      be1=-pcm(1)/pcm(4)
      be2=-pcm(2)/pcm(4)
      be3=-pcm(3)/pcm(4)
      ss=pcm(4)**2-(pcm(1)**2+pcm(2)**2+pcm(3)**2)
      if(ss.le.0.0d0) goto 110
      gam=pcm(4)/sqrt(ss)
 
      do i=1,4
       cu1(i)=cu(i)
       cub1(i)=cub(i)
      do j=1,4
        ten1(i,j)=ten(i,j)
      end do
      end do
 
      call jamrobo(0.0d0,0.0d0,be1,be2,be3,gam,cu1(1),cu1(2),cu1(3), 
     & cu1(4))
      call jamrobo(0.0d0,0.0d0,be1,be2,be3,gam
     $                               ,cub1(1),cub1(2),cub1(3),cub1(4))
      call jamlemt(be1,be2,be3,gam,ten1)
 
      a(ie,1)=cu1(4)                   ! density
      a(ie,2)=cub1(4)                  ! baryon density
      a(ie,3)=ten1(4,4)                ! energy density
      a(ie,4)=(ten1(1,1)+ten1(2,2))/2  ! (t_xx+t_yy)/2
      a(ie,5)=ten1(3,3)                ! t_zz

 110  continue
c...Local frame 2
      bex=-cu(1)/cu(4)
      bey=-cu(2)/cu(4)
      bez=-cu(3)/cu(4)
      gamm=cu(4)/sqrt(cc)
      call jamrobo(0.0d0,0.0d0,bex,bey,bez,gamm,cu(1),cu(2),cu(3),cu(4))
      call jamrobo(0.0d0,0.0d0,bex,bey,bez,gamm,cub(1),cub(2),cub(3), 
     & cub(4))
      call jamlemt(bex,bey,bez,gamm,ten)

      a(ie,6)=cu(4)                  ! density
      a(ie,7)=cub(4)                 ! baryon density
      a(ie,8)=ten(4,4)               ! energy density
      a(ie,9)=(ten(1,1)+ten(2,2))/2  ! (t_xx+t_yy)/2
      a(ie,10)=ten(3,3)              ! t_zz

c-----------------------------------------------------------------------
 
c...a(1)  density
c...a(2)  baryon density
c...a(3)  energy density
c...a(4)  (t_xx+t_yy)/2
c...a(5)  t_zz

c...a(6)  density
c...a(7)  baryon density
c...a(8)  energy density
c...a(9)  (t_xx+t_yy)/2
c...a(10) t_zz

c...a(11) Lorentz invariant density
c...a(12) Lorentz invariant baryon density
c...a(13) Lorentz invariant energy density
c...a(14) Lorentz invariant pressure

c-----------------------------------------------------------------------

          tempt1(ie)=0.0d0
          tempz1(ie)=0.0d0
          if(a(ie,1).gt.0.0d0) then
            tempt1(ie)=a(ie,4)/a(ie,1)
            tempz1(ie)=a(ie,5)/a(ie,1)
          endif

          tempt2(ie)=0.0d0
          tempz2(ie)=0.0d0
          if(a(ie,6).gt.0.0d0) then
            tempt2(ie)=a(ie,9)/a(ie,6)
            tempz2(ie)=a(ie,10)/a(ie,6)
          endif

          temp(ie)=0.0d0
          if(a(ie,11).gt.0.0d0) then
            temp(ie)=a(ie,14)/a(ie,11)
          endif

          a(ie,1)=a(ie,1)/rho0
          a(ie,2)=a(ie,2)/rho0
          a(ie,6)=a(ie,6)/rho0
          a(ie,7)=a(ie,7)/rho0
          a(ie,11)=a(ie,11)/rho0
          a(ie,12)=a(ie,12)/rho0

 200      continue
          write(ida,810)it*parc(7)
     $    ,(a(ie,j),j=1,3),tempt1(ie),tempz1(ie)
     $    ,(a(ie,j),j=6,8),tempt2(ie),tempz2(ie)
     $    ,(a(ie,j),j=11,14),temp(ie)

          end do  ! end time loop)
            call jamfile(1,166,' ')
          end do  ! end ie loop)


      endif

800   format('# Particle/energy density box=1*1*1'
     $  /'# time  rho rho_b e(GeV/fm^3) T_xy  T_z')
801   format('# Particle/energy density  vol=2piR*z'
     $  /'# time  rho rho_b e(GeV/fm^3) T_xy  T_z')
810   format(f8.4,1x,3(f9.3,1x),2(f9.3,1x)
     $              ,3(f9.3,1x),2(f9.3,1x)
     $              ,4(f10.3,1x),f10.3)
 
      end

c***********************************************************************

      subroutine jamadns2(msel)

c...Purpose: calculation of density using Gauss smearing
      include 'jam1.inc'
      include 'jam2.inc'
       
      parameter(maxt=100,mab=9)
      parameter(rho0=0.168d0)
      parameter(dr=0.5d0,dh=0.2d0,dz=0.2d0,zz=1.0d0)
c     parameter(widg=2*0.6d0)                   ! Gaussian width
      parameter(widg=2*1.0d0)                   ! Gaussian width

      dimension a(2,mab)
      dimension g(4,4),cu(4),cub(4),ten(4,4)
      dimension cur(2,4,maxt),curb(2,4,maxt),tens(2,4,4,maxt)
      dimension tempt(2),tempz(2),temp(2)
      dimension a2(mab),b2(mab),a3(mab)
      save kt,cur,curb,tens,ktmax
      data g/ -1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,1/
      data istime/1/
c=======================================================================

c...Reset counters
c-----------------------------------------------------------------------
      if(msel.eq.0) then

      ktmax=0
      do it = 1,maxt

      do i=1,4
        cur(1,i,it)=0.0d0
        cur(2,i,it)=0.0d0
        curb(1,i,it)=0.0d0
        curb(2,i,it)=0.0d0
      do j=1,4
        tens(1,i,j,it)=0.0d0
        tens(2,i,j,it)=0.0d0
      end do
      end do

      end do

c-----------------------------------------------------------------------
      else if(msel.eq.1) then

        kt=0

c-----------------------------------------------------------------------
      else if(msel.eq.2) then

c     b=pard(2)
c     r1=pard(40)
c     r2=pard(50)
c     rx1=min(r1+b/2,r2-b/2)
c     rx2=min(r1-b/2,r2+b/2)
c     rr2=(min(rx1,rx2))**2
      rr2=(min(pard(40),pard(50)))**2

c     maxz=nint(zz/dz)
c     maxr=nint(sqrt(rr2)/dr)
c     maxh=nint(1.0/dh)
c     maxm=(maxz+1)*(maxr+1)*(maxh+1)
      maxz=0
      maxr=0
      maxh=0
      maxm=1

      time=pard(1)
      if(istime.eq.1) then
        kt=kt+1
        if(kt.gt.maxt) return
        ktmax=max(kt,ktmax)
      endif

      widcof=(1.d0/(3.14159d0*widg))**1.5d0

c==============================================================
c     do ir=0,maxr
c     do ih=0,maxh
c     do iz=0,maxz
c==============================================================

c     x2=ir*dr*cos(ih*dh)
c     y2=ir*dr*sin(ih*dh)
c     z2=iz*dz
      x2=0.0d0
      y2=0.0d0
      z2=0.0d0

c---------------------------------------------------------------------
c....Loop over all particles
      do 100 i=1,nv
 
        k1=k(1,i)
        if(k1.gt.10) goto 100   ! dead particle
        if(p(5,i).le.1d-5) goto 100
c       if(r(5,i).gt.time) goto 100
        if(abs(k(7,i)).eq.1) goto 100   ! not yet interaction

        dt=time-r(4,i)
        if(dt.lt.0.0d0) goto 100

        x1=r(1,i)+dt*p(1,i)/p(4,i)-x2
        y1=r(2,i)+dt*p(2,i)/p(4,i)-y2
        z1=r(3,i)+dt*p(3,i)/p(4,i)-z2

c...Proper time
        if(istime.eq.2) then
          tau2=time**2-z1**2
          if(tau2.le.0.0d0) goto 100
          tau=sqrt(tau2)
          kt=int(tau/parc(7))
          if(kt.le.0.or.kt.gt.maxt) goto 100
          ktmax=max(kt,ktmax)
        endif

c       bar=kfprop(k(2,i),2)/3.0d0
        bar=k(9,i)/3.0d0
        if(r(5,i).gt.time) then
          iq=mod(abs(k(1,i))/10,10)
          if(iq.eq.3) iq=2
          bar=iq/3d0*isign(1,k(2,i))
          if(k(9,i).eq.0)bar=abs(bar)
        endif

        xtra=x1**2+y1**2+z1**2
     $          +((x1*p(1,i)+y1*p(2,i)+z1*p(3,i))/p(5,i))**2
        if(xtra/widg.gt.30.d0) goto 100

c       nc1=nc1+1
        gam=p(4,i)/p(5,i)
        den=widcof*gam*exp(-xtra/widg)

        do im=1,4
          cur(1,im,kt)=cur(1,im,kt) + p(im,i)/p(4,i)*den
          curb(1,im,kt)=curb(1,im,kt) + p(im,i)/p(4,i)*den*bar
        do in=1,4
          tens(1,im,in,kt)=tens(1,im,in,kt) + p(im,i)*p(in,i)/p(4,i)*den
        end do
        end do

        if(r(5,i).gt.time) goto 100

c       nc2=nc2+1
        do im=1,4
          cur(2,im,kt)=cur(2,im,kt) + p(im,i)/p(4,i)*den
          curb(2,im,kt)=curb(2,im,kt) + p(im,i)/p(4,i)*den*bar
        do in=1,4
          tens(2,im,in,kt)=tens(2,im,in,kt) + p(im,i)*p(in,i)/p(4,i)*den
        end do
        end do

100   continue

c...Output
      else if(msel.eq.3) then
c---------------------------------------------------------------------
        wei=1.0d0/dble(mstd(21)*mstc(5))
        ida=mstc(36)
        call jamfile(0,167,' ')
        write(ida,800)mstd(21),widg/2
        write(ida,805)
c...Loop over time.
      do 200 it=1,ktmax

        do j=1,mab
        a(1,j)=0.0d0
        a(2,j)=0.0d0
        end do

      do 250 in=1,2

c       if(in.eq.1.and.nc1.eq.0) goto 250
c       if(in.eq.2.and.nc2.eq.0) goto 250


c     do i=1,mab
c      a2(i)=0.0
c      b2(i)=0.0
c     end do

      do i=1,mab
       a3(i)=0.0d0
      end do

      do i=1,4
        cu(i)=cur(in,i,it)*wei
        cub(i)=curb(in,i,it)*wei
      do j=1,4
        ten(i,j)=tens(in,i,j,it)*wei
      end do
      end do

c...Lorentz invariant Density
      cc=cu(4)**2-(cu(1)**2+cu(2)**2+cu(3)**2)
      if(cc.gt.0.0d0) then
        a3(6)=sqrt(cc)
      else
        goto 250
      endif

c...Lorentz invariant Baryon density
      bnorm=cub(4)**2 -( cub(1)**2+cub(2)**2+cub(3)**2 )
      if(bnorm.gt.0.0d0) then
        a3(7)=sqrt(bnorm)
      endif

c...Lorentz invariant pressure and energy density
      a3(8)=0.0d0
      a3(9)=0.0d0
      do i=1,4
      do j=1,4
       tmp=g(i,j)*cu(i)*cu(j)/cc
       a3(8)=a3(8)+ten(i,j)*tmp                    ! energy density
       a3(9)=a3(9)-1.d0/3.d0*ten(i,j)*(g(i,j)-tmp)     ! pressure
      end do
      end do

c     a3(1)=cu(4)                  ! density
c     a3(2)=cub(4)                 ! baryon density
c     a3(3)=ten(4,4)               ! energy density
c     a3(4)=(ten(1,1)+ten(2,2))/2  ! t_xx
c     a3(5)=ten(3,3)               ! t_zz

c...Local frame
      be1=-cu(1)/cu(4)
      be2=-cu(2)/cu(4)
      be3=-cu(3)/cu(4)
      gam=cu(4)/sqrt(cc)
      call jamrobo(0.0d0,0.0d0,be1,be2,be3,gam,cu(1),cu(2),cu(3),cu(4))
      call jamrobo(0.0d0,0.0d0,be1,be2,be3,gam,cub(1),cub(2),cub(3), 
     & cub(4))
      call jamlemt(be1,be2,be3,gam,ten)
      a3(1)=cu(4)                  ! density
      a3(2)=cub(4)                 ! baryon density
      a3(3)=ten(4,4)               ! energy density
      a3(4)=(ten(1,1)+ten(2,2))/2  ! t_xx
      a3(5)=ten(3,3)               ! t_zz

c...Max. and average
c     if(in.eq.1) then
c     do j=1,mab
c       a2(j)=max(a2(j),a3(j))
c     end do
c     else
c     do j=1,mab
c       b2(j)=max(a2(j),a3(j))
c     end do
c     endif

c================================================================
c     end do
c     end do
c     end do
c================================================================

      do j=1,mab
        a(in,j)=a(in,j)+a3(j) 
      end do

 250  continue

c-----------------------------------------------------------------------

c...a(1)  density
c...a(2)  baryon density
c...a(3)  energy density
c...a(4)  (t_xx+t_yy)/2
c...a(5)  t_zz
c...a(6)  Lorentz invariant density
c...a(7)  Lorentz invariant baryon density
c...a(8)  Lorentz invariant energy density
c...a(9)  Lorentz invariant pressure

c-----------------------------------------------------------------------

          do ie=1,2

          tempt(ie)=0.0d0
          tempz(ie)=0.0d0
          if(a(ie,1).gt.0.0d0) then
            tempt(ie)=a(ie,4)/a(ie,1)
            tempz(ie)=a(ie,5)/a(ie,1)
          endif

          temp(ie)=0.0d0
          if(a(ie,6).gt.0.0d0) then
            temp(ie)=a(ie,9)/a(ie,6)
          endif

          a(ie,1)=a(ie,1)/rho0
          a(ie,2)=a(ie,2)/rho0
          a(ie,6)=a(ie,6)/rho0
          a(ie,7)=a(ie,7)/rho0

          end do

          write(ida,810)it*parc(7)
     $ ,(a(1,j),j=1,3),tempt(1),tempz(1),(a(1,j),j=6,9),temp(1)
     $ ,(a(2,j),j=1,3),tempt(2),tempz(2),(a(2,j),j=6,9),temp(2)

 200  continue

        call jamfile(1,167,' ')

      endif

800   format('# Particle/energy density (Gauss smear)'
     $  /'# Event=',i5,' Gaussian width=',f8.3
     $  /'# (1) time'
     $  /'# (2) rho                         (12)'
     $  /'# (3) rho_b                       (13)'
     $  /'# (4)  e(GeV/fm^3)                (14)'
     $  /'# (5)  T_xy                       (15)'
     $  /'# (6)  T_z                        (16)'
     $  /'# (7)  Lorentz invariant rho      (17)'
     $  /'# (8)  Lorentz invariant rho_b    (18)'
     $  /'# (9)  Lorentz invariant e        (19)'
     $  /'# (10) Lorentz invariant p        (20)'
     $  /'# (11) Lorentz invariant T        (21)')

805   format('#',3x,'(1)',7x,'(2)',7x,'(3)',7x,'(4)',8x,'(5)',6x,'(6)'
     $  ,8x,'(7)',8x,'(8)',8x,'(9)',8x,'(10)',7x,'(11)'
     $  ,6x,'(12)',6x,'(13)',6x,'(14)',6x,'(15)',7x,'(16)',5x,'(17)'
     $  ,7x,'(18)',7x,'(19)',7x,'(20)',7x,'(21)',7x,'(21)')
810   format(f8.4,1x,3(f9.3,1x),2(f9.3,1x),4(f10.3,1x),f10.3
     $           ,1x,3(f9.3,1x),2(f9.3,1x),4(f10.3,1x),f10.3)
 
      end

c***********************************************************************

      subroutine jamlemt(be1,be2,be3,gamma,emt)

c...Lorentz transformation of energy-momentum tensor.

      implicit double precision(a-h, o-z)
      dimension emt(4,4),a(4,4),emtr(4,4),beta(3)

      beta(1)=-be1
      beta(2)=-be2
      beta(3)=-be3
      eta=gamma**2/(1+gamma)

      a(4,4)=gamma
      do j=1,3
        a(j,4)=-gamma*beta(j)
        a(4,j)=-gamma*beta(j)
        do i=1,3
          if (i.eq.j) then
            a(i,j) = 1+eta*beta(i)**2
          else
            a(i,j) = eta*beta(i)*beta(j)
          endif
        enddo
      enddo

      do j=1,4
      do i=1,4
        emtr(i,j)=0 
      enddo
      enddo

      do j=1,4
      do i=1,4
      do n=1,4
      do m=1,4
        emtr(i,j)=emtr(i,j)+a(i,m)*a(j,n)*emt(m,n)
      enddo
      enddo
      enddo
      enddo
 
      do j=1,4
        do i=1,4
          emt(i,j)=emtr(i,j)
        enddo
      enddo

      end

c***********************************************************************

      subroutine jameksve(msel)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Remember all e_kin of particles which will be used for temperature.
c   iset: specify the memory set e.g. time step.
c   nrem: number of particles memorized in the set.
c   ek:   memorized kinetic energy of particles.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      include 'jam1.inc'
      include 'jam2.inc'

      parameter (maxset=100,maxrem=2000,maxnbin=50)
      common/edist/ek(maxset,maxrem),count(maxset,0:maxnbin),de(maxset),
     &             s_time2(maxset),nrem(maxset),nset,nbin(maxset)
      save /edist/
      save iset0


      if(msel.eq.1) then
        iset0=0
        return
      endif

c...Remember ek 

      iset0=iset0+1

      if(iset0.gt.nset)nset=iset0
      s_time2(iset0)=pard(1)

c...Remember only particles alive
      nn=0
      do i=1,nv
        if(k(1,i).le.11) then
          if(nrem(iset0).lt.maxrem) then
            nn=nn+1
            ekin=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)-p(5,i)
            nrem(iset0)=nrem(iset0)+1
            ek(iset0,nrem(iset0))=ekin
          else
            call jamerrm(1,0,'(jameksve:)number to remember exceeds')
          endif
        endif
      enddo

      write(6,*) 'rem',iset0,'/',nset ,'  nrem=',nn,'/',nrem(iset0)
      return

c***********************************************************************

      entry jamekclr(iset1)

c...Clear ek memory.

      if(iset1.eq.0)then
        do iset=1,maxset
          do i=1,maxrem
            ek(iset,i)=0.d0
          enddo
          nrem(iset)=0
        enddo
      else
        iset=iset1
        do i=1,maxrem
          ek(iset,i)=0.d0
        enddo
        nrem(iset)=0
        nset=0
      endif
      return
      end

c***********************************************************************

      subroutine jamtout

c...Output temperature in box calculation.
      include 'jam1.inc'
      include 'jam2.inc'
      parameter (maxset=100,maxrem=2000,maxnbin=50)
      common/edist/ek(maxset,maxrem),count(maxset,0:maxnbin),de(maxset),
     &             s_time2(maxset),nrem(maxset),nset,nbin(maxset)
      save /edist/
      dimension temp(maxset)

      idp=mstc(36)

c...Print temperature by time.
      call jamfile(0,168,' ')
      
        write(idp,100)
1000    format('Xmin 0;Ymin 0',/,'def time(1) t% 1 1 1 0 1 1 &1 fm/c')
1100    format('c;time ',f5.0)
100     format('#time(fm/c) temperature N')
200     format(f8.4,1x,d10.4,1x,i10)
300     format('# ',f8.4,1x,e10.4,1x,i10,/,'Lh 1;c')
400     format(f8.4,1x,e10.4)
500     format('L 1;c 12',/,'rep ee ',3(f8.4,1x),
     &         ' ee @temp=',e9.4,',2*sqrt(ee/pi/temp^3)*exp(-ee/temp)')
600     format('a 10 ',2f5.0) 
610     format('a% 0 -1 1 ') 
620     format('a% 1.1 ',i3,' 1 ') 

      do iset=1,nset
        call jamctmp(iset,temp(iset))
        write(idp,200)s_time2(iset),temp(iset),nrem(iset)
      enddo
      call jamfile(1,168,' ')

c...Print e_kin distribution (temperature fit)
      if(mstc(169).eq.1) then
        ratio=2.7d0 ! graph distribution parameter
        nord=int(sqrt(nset/ratio)+0.5d0)
        nvert=int(1.0d0*nset/nord+0.999d0)
        height=240/nvert
        call jamfile(0,169,' ')
        write(idp,1000)
        do iset=1,nset
          if(iset.eq.1) then
            write(idp,600) 260-height,height
          else if(mod(iset-1,nvert).eq.0) then
            write(idp,620) nvert-1
          else
            write(idp,610)
          endif
          write(idp,300)s_time2(iset),temp(iset),nrem(iset)
          do ie=1,nbin(iset)
            write(idp,400)ie*de(iset),count(iset,ie)
          enddo
          write(idp,500)de(iset),nbin(iset)*de(iset),de(iset),temp(iset)
          write(idp,1100) s_time2(iset)
        enddo
        call jamfile(1,169,' ')
      endif

      end

c***********************************************************************

      subroutine jamctmp(iset,temperature)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   iset: specify the memory set e.g. time step.
c   nrem: number of particles memorized in the set.
c   ek:   memorized kinetic energy of particles.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      include 'jam1.inc'
      include 'jam2.inc'

      parameter (maxset=100,maxrem=2000,maxnbin=50)
      common/edist/ek(maxset,maxrem),count(maxset,0:maxnbin),de(maxset),
     &             s_time2(maxset),nrem(maxset),nset,nbin(maxset)
      save /edist/


c...Calculate temperature from memorized ek

c...Determine the enery bin so that empty bin does not exist.
      edifmax=0.d0
      ekinmax=ek(iset,1)
      ekinmin=ek(iset,1)
      do i=2,nrem(iset)
c....Determine e_kin max and min.
        ekin=ek(iset,i)
        if(ekinmax.lt.ekin) then
          ekinmax=ekin
          imaxkin=i
        endif
        if(ekinmin.gt.ekin) then
          ekinmin=ekin
          iminkin=i
        endif
c...Determine minimum diff of ek for each i.
        edifmin=ekin
        do j=1,nrem(iset)
          ekin2=ek(iset,j)
          if(ekin.gt.ekin2.and.ekin-ekin2.le.edifmin) then
            edifmin=ekin-ekin2
            jmin=j
          endif
        enddo
c.....determine max of minimum diff of ek.
        if(edifmin.gt.edifmax) then
          edifmax=edifmin
        endif
      enddo

c...first way of bin.
      if(mstc(168).eq.2) then
        nbin(iset)=int(ekinmax/edifmax)
        if(nbin(iset).gt.maxnbin)nbin(iset)=maxnbin
c...another way of bin.
      else if(mstc(168).eq.1) then
        nbin(iset)=maxnbin
      endif


c....now bin is determined.

      de(iset)=ekinmax/nbin(iset)
      do i=0,nbin(iset)
        count(iset,i)=0.d0
      enddo

      do i=1,nrem(iset)
        ekin=ek(iset,i)
        ibin=nint(ekin/de(iset))
        if(ibin.le.nbin(iset)) count(iset,ibin)=count(iset,ibin)+1.d0
      enddo
      do i=0,nbin(iset)
        if(mstc(168).eq.2) then
          if(count(iset,i).lt.0.5d0) then
            write(check(1),'(3(i9,1x))')i,nbin(iset),count(iset,i)
            call jamerrm(1,1,'some count has empty value.')
          endif
        endif
        count(iset,i)=count(iset,i)/nrem(iset)/de(iset)
      enddo
c     do i=0,nbin(iset)
c      write(*,*)i*de(iset),count(iset,i)
c     enddo

      ntry=5000
      smin=1000000000.d0
      pi=3.141593d0
      do itry=1,ntry
        temp=pard(22)*2*itry/ntry
        s0=0
        do i=1,nbin(iset)
          ee=i*de(iset)
          s0=s0+(2*sqrt(ee/pi/temp**3)*exp(-ee/temp)-count(iset,i))**2
        enddo
        if(s0/nbin(iset).lt.smin) then
          smin=s0/nbin(iset)
          tempmin=temp
        endif
      enddo
      temperature=tempmin
      chi2=smin
      end
      


c***********************************************************************

      subroutine jamgout(msel,rr,in)

c...Booking of ground state nucleus.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      parameter (nmax=50,dr=0.2d0)
      dimension a(2,nmax)

c...Initialize
      if(msel.eq.0) then
         do i=1,nmax
           a(2,i)=0
         end do

c...Booking
      else if(msel.eq.2) then

         ix=rr/dr
         if(ix.le.0.or.ix.gt.nmax) return
         a(in,ix)=a(in,ix)+1
       
c...Output
      else if(msel.eq.3) then

c....Open file
        idb=mstc(36)
        call jamfile(0,157,' ')
        write(idb,800)
800   format('# Ground state of target nucleus: r  targ. proj.')

        do ix=1,nmax
          xx=dr*(ix-0.5d0)
          write(idb,810)xx,(a(j,ix)/dble(mstd(21)),j=1,2)
        end do
810     format(f9.4,6(1x,f9.4))

c...Close file
        call jamfile(1,157,' ')

      endif

      end
