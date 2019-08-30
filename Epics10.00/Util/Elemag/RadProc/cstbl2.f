!     ****************************************************************
!     *                                                              *
!     * cstbl$2:   create sampling table for brems and pair          *
!     *            with landau effect at high energies               *
!     *                                                              *
!     ****************************************************************
!
!      prepare:   fbrem in fbrem2 module ( -inc works for this end)
!
!                 others same as cstbl$1
!
!
!        sub0  fbrem,         smigb, sbrem2,gzai, fbremi, bremr,mkdt,
!     rutbnr, gquadt,gmigdl, psimig, fbremc, gqdir1
!
!     lines between  ---  should be used in sampling routines too
!
!      -------------------- landau ------------------------------
!          **** electron ****
!
      data egmnh/1.e-4/
      parameter (eeh1=20., eeh2 = 5000.e3, neeh=25, dueh1=.1,
     *  ueh1=0., ueh2=.7, nueh1=(ueh2-ueh1+.00001)/dueh1+1)
      dimension  tcbh(neeh)
      dimension  be1(nueh1, neeh)
!
!                         sqrt(1-.7)
      parameter (ueh3=0., ueh4=0.5477225,
     *nueh2=11, dvueh2=nueh2-1, dueh2=(ueh4 - ueh3)/dvueh2)
      dimension be2(nueh2,neeh)
!
!
!         *** gamma ****
!
       parameter ( egh1=10000., egh2= 5000.e3, negh=20, dugh2=.1,
     * ugh3=0., ugh4=1., nugh2=(ugh4-ugh3+.00001)/dugh2+1)
      dimension tcph(negh)
      dimension vp(nugh2, negh)
!
!
!     ---------------------------------------------------------
!
      dimension eel(neeh), eell(neeh)
      dimension egm(negh), egml(negh)
!
      character matter*4, cap*100, media*45, mem*6, file*25
!
      data  eps/0.001/
      common /cbrem/upsi,emain
      external bremr,pairr
!
      common /landuc/ x0cm, x0g, s1, alogs1, sconst
      common /landu1/e
!
      dimension ua(100)
!
!
!
!      material   z     a     rho     t0 (g/cm2)  t0 (cm)    ec
!         pb      82   207.2  11.35    5.82        0.513     6.7
!         cu      29    63.5   8.94   12.7        1.42     16+x
!         fe      26    55.85  7.86   13.7         1.74     18+x
!          w      74   183.92  19.3    6.275      .325
!
      jpunch=0
      io=0
      write(6,'('' enter matter'')')
      read(5, *) matter
      if(matter .eq. 'pb') then
          z=82.
          a=207.2
          rho=11.35
      elseif(matter .eq. 'fe') then
          z=26.
          a=55.85
          rho=7.86
      elseif(matter .eq. 'w') then
          z=74.
          a=183.92
          rho=19.3
      elseif(matter .eq. 'cu') then
          z=29.
          a=63.5
          rho=8.94
      else
          write(6, '('' matter='',a4, '' undefined'')')
          write(6, '('' enter matter, z, a, rho '')')
          read(5, *) matter, z, a, rho
      endif
!
      write(6,'('' media='',a4, '' z='', f6.2, '' a='',f6.2,
     * '' rho='',g13.3 )') matter, z, a, rho
!      write(media,'('' media='',a4, '' z='', f6.2, '' a='',f6.2,
!     * '' rho='',g9.3 )') matter, z, a, rho
!
      call zpart( z, a, rho)
      write(6,'('' gamma cut off energy ='', g9.3,''gev'')') egmnh
!
      lmx=6
      jpunch = 0
!
      deelh=log10(eeh2/eeh1)/(neeh-1)
      deghl=log10(egh2/egh1)/(negh-1)
      
       do   l=6, lmx
          write(l, '(''      data deelh/'',f9.6,''/'')') deelh
          write(l, '(''      data deghl/'',f9.6,''/'')') deghl
       enddo
!
!
      de=10.**deelh
      e=eeh1
       do   ie=1,neeh
          eel(ie)=e
          eell(ie)=log10(e)
          vmin=egmnh/e
          call totcb(vmin,1.,ans)
          tcbh(ie)=ans
          emain=e
          vmax=1.
          tp=tcbh(ie)
          u=ueh1
          vl = vmin
          do   iu=2,nueh1
              u=u+dueh1
              upsi=u*tp
              call kbchop(bremr, vl, vmax, eps, v, j)
              if(j .le. 0) then
                 write(*, *) ' error; e=',e,' u =',u, ' v =',v
              endif
              be1(iu,ie)= log(v/vmin)
           enddo
          be1(1 ,ie)=   log(1./vmin)
          e=e*de
       enddo
!
!
!
      u=ueh1
       do   i=1, nueh1
         ua(i)=u
         u=u+dueh1
       enddo
!
!
!
       do   l=6,lmx
           el1=alog10(eeh1)
           el2=alog10(eeh2)
           write(l, '(''c'')')
           write(l,'(''c          energy of electron''/
     *     (1hc,5x, 5g12.4))') eel
           write(l,'(''c         log10 step='',f9.6)') deelh
           write(l, '(''c'')')
           write(l,'(''c         log10 of electron energy boundary'')')
           write(l,'(''      data eeh1l/'',f9.6,''/'')') el1
           write(l,'(''c    *     ,eeh2l/'',f9.6,''/'')') el2
           write(l, '(''c'')')
       enddo
       call kmkDataStm1(tcbh, neeh, 'tcbh', 'f8.4', 8)
!        call mkdt('tcbh', tcbh, 1, neeh,      'f8.4,   ', 8,jpunch)
!
       do   l=6, lmx
           write(l, '(''c'')')
           write(l,'(''c          be1(iu,ie)= alog(v/vmin)'')')
           write(l,'(''c          from u='',f5.3,'' to '',f5.3,
     *     '' step='',f5.3)') ueh1, ueh2, dueh1
           write(l,'(''c          from log10(e)='',f8.4,'' to'',f8.4,
     *     '' step'', f8.4)') el1, el2, deelh
           write(l,'(''c          dim. of u='',i3, '' dim of e='',i3)')
     *     nueh1,neeh
           write(l, '(''c'')')
       enddo
       call kmkDataStm1(be1, nueh1*neeh, '  be1', 'f9.5', 9)
!        call mkdt('   be1', be1,   1, nueh1*neeh, 'f9.5,    ',9,jpunch)
!
      e=eeh1
      do   ie=1,neeh
          vmin=egmnh/e
          sqrtv=sqrt(vmin)
          emain=e
          vl=vmin
          vmax=1.
          tp= tcbh(ie)
          us=ueh4
           do   iu=nueh2,2,-1
              u=1. -us**2
              upsi=u*tp
              call kbchop(bremr, vl, vmax, eps, v, j)
              if(j .le. 0) then
                 write(*, *) ' error for be2, e=',e, ' u=',u
              endif
              be2(iu,ie)= (sqrt(v) -sqrtv)/(1.-u)
              us=us-dueh2
           enddo
!
          be2(1,ie)=tp /2/fbrem(vmin)/sqrtv
          e=e*de
       enddo
!
!
       us=ueh3
       do   i=1, nueh2
         ua(i)=us
         us=us+dueh2
       enddo
!
!

       do   l=6,lmx
           write(l,'(''c'')')
           write(l, '(''c         be2(iu,ie)= (sqrt(v) -sqrtv)/(1.-u)''
     *     )')
           write(l,'(''c          from sqrt(1-u)='',f5.3,'' to '',f5.3,
     *     '' step='', f6.3)') ueh3, ueh4, dueh2
           write(l,'(''c          from log10(e)='',f9.6,'' to'',f9.6,
     *     '' step'', f8.4)') el1, el2, deelh
           write(l,'(''c'')')
       enddo
       call kmkDataStm1(be2, nueh2*neeh, '  be2', 'f8.4', 8)
!c      call mkdt('   be2',be2, 1, nueh2*neeh, 'f8.4,    ',8,jpunch)
!
!
!
!          gamma
!
      e=egh1
      de=10.**deghl
!
       do   ie=1,negh
          egm(ie)=e
          egml(ie)=log10(e)
          call totcp(.5, 1., ans)
          tcph(ie)=ans*2
          tp=ans
          emain=e
          u=ugh3
          vmax=1.
          vl=.5
          do   iu=2,nugh2-1
              u=u+dugh2
              upsi=u*tp
              call kbchop(pairr, vl, vmax, eps, v, j)
              if(j .le. 0) then
                 write(*,*) ' error for vp, e=',e,' u=',u
              endif
              vp(iu,ie)=v
           enddo
          vp(1,ie)=1.
          vp(nugh2,ie)=0.5
          e=e*de
       enddo
!
!
!
      u=ugh3
       do   i=1, nugh2
         ua(i)=u
         u=u+dugh2
       enddo
!
!
!
!
       do   l=6, lmx
           el1=alog10(egh1)
           el2=alog10(egh2)
           write(l,'(''c'')')
           write(l,'(''c'')')
           write(l,'(''c         energy of gamma''/
     *     (1hc,5x, 5g12.4))') egm
           write(l,'(''c       log10 step='',f9.6)') deghl
           write(l,'(''c'')')
           write(l,'(''c       log10 of gamma energy boundary'')')
           write(l,'(''      data egh1l/'',f9.6,''/'')') el1
           write(l,'(''c    *     ,egh1l/'',f9.6,''/'')') el2
           write(l,'(''c'')')
       enddo
       call kmkDataStm1(tcph, negh, ' tcph', 'f8.4', 8)
!      call mkdt('tcph', tcph, 1, negh,      'f8.4,   ',8, jpunch)
!
!
       do   l=6, lmx
           write(l,'(''c'')')
           write(l,'(''c        table of v'')')
           write(l,'(''c        from u='',f7.4,'' to '',f7.4,'' step'',
     *     f8.4)') ugh3, ugh4, dugh2
           write(l,'(''c        from log10(e)='',f9.6,'' to'', f9.6,
     *     '' step'', f8.4)') el1, el2, deghl
           write(l,'(''c'')')
       enddo
       call kmkDataStm1(vp, negh*nugh2, '  vp', 'f8.4', 8)
!      call mkdt('   vp',  vp,    1,  negh*nugh2, 0.,    0, jpunch)
      end
!     ****************************************************************
!     *                                                              *
!     * bremr:  used to solve equation for making sampling table     *
!     *         for brems                                            *
!     *                                                              *
!     * pairr: //  for pair                                          *
!     *                                                              *
!     ****************************************************************
!
!
      function bremr(v)
!
      common /cbrem/upsi,emain
      call totcb(v,1.,ans)
      bremr=ans/upsi-1.
      return
!
!     ***********
      entry pairr(v)
!     ***********
!
      call totcp(v , 1., ans)
      bremr=ans/upsi - 1.
      return
      end
!     *****************
      subroutine totcb(vmin,vmax,ans)
!     *****************
!
!        integration of bremsung function from vmin to =max.
!
      external fbrem,fpair
!
      ans=0.
      v2=vmax
    5 continue
      v1=v2/10.
      v1=amax1(v1,vmin)
      call gquadt(fbrem,v1,v2,ans1)
      ans=ans+ans1
      if(v1 .eq. vmin) return
      v2=v1
      goto 5
!
!
!
!     ***********
      entry totcp(vmin,vmax,ans)
!     ***********
!
!
!        integralation of pair-cre function
!
!
!
      d=(vmax-vmin)/30.
      vt=vmax-d
      call gquadt(fpair,vmin,vt,ans1)
      call gquadt(fpair,vt,vmax,ans2)
      ans=ans1+ans2
      end
!     *************************************  not used ****************
!     *                                                              *
!     * fbremc:  sqrt(v)*fbrem(v):                                   *
!     *            with landau effect at high energies               *
!     *                                                              *
!     ****************************************************************
!
!
!
      function fbremc(v)
!
!
!
!        sqrt(v)*fbrem(v):   at v=0, equals to 8*sqrt(2 * sconst/e)
!        =264.3/sqrt(e) for pb
!
      common/landu1/e
      common/landu2/vcut,f0
!
      common /landuc/ x0cm, x0g, s1, alogs1, sconst
!
    5 continue
      if(v .ne. 0.) goto 10
      fbremc=8. *  sqrt( 2. * sconst /e )
      return
   10 continue
      fbremc=sqrt(v)*fbrem(v)
      return
!
!
!     ************
      entry fbrems( v )
!     ************
!
      if(v .ge. vcut) goto 5
      fbremc=f0
      return
      end
!     ********************************************  not used  ********
!     *                                                              *
!     * fbremi:  integral of brems function                          *
!     * fpairi:  //          pair                                    *
!     * tprbb:   total probability of brems                          *
!     *                                                              *
!     ****************************************************************
!
!
!        integral of fbrem from v to 1. at energy ee(gev).
!        for a given ee, you must use tprbb before using fbremi.  fpair
!        can be called at any time.  total prob may be given by fpairi(
!        e,.5)*2.
!
!
!
      function fbremi(ee,v)
!
!
!
      common/landu1/e
      data v0,v1,v2/ 0.008, 0.05, 0.1 /
      data e1, e2/ 3000., 100. /
      external fbremc,fpair
      e=ee
      if( e .gt. e1) goto 10
      if(e .gt. e2) goto 20
      ans=cibf(v0,v2,v)
      goto 50
   10 continue
      call gqdir2( fbremc, 0., v, 2, ans)
      goto 50
   20 continue
      ans=cibf(v1,1.,v)
   50 continue
      fbremi=tprbbs-ans
      return
!
!
      entry fpairi(ee,v)
!
!
      e=ee
      call gquadt(fpair,v,1.,ans)
      fbremi=ans
      return
!
!
      entry tprbb(ee)
!
!        total probability of bremsung at energy ee(gev)
      e=ee
      if(e .gt. e1) goto 110
      if(e .gt. e2) goto 120
      tprbbs=cibf(v0,v2,1.)
      goto 150
  110 continue
      call gqdir2(fbremc, 0., 1., 2, tprbbs)
      goto 150
  120 continue
      tprbbs=cibf(v1,1.,1.)
  150 continue
      fbremi=tprbbs
      return
      end
!     ****************************************************  not used *
!     *                                                              *
!     * cibf:  auxliary function to integrate brems function         *
!     *        with landau effect                                    *
!     *                                                              *
!     ****************************************************************
!
!        integral of bremsung function when landau effect exists.  inte
!        gral of bremsung function (f) from 0 to v is given by using gq
!        ir2.   if v is .gt. v0, or v1, integral is divided 2 or 3 port
!        on. ( at low wenery f is steep at small v, so that it is neces
!        ary to divide integrals.  energy must be informed thru common/
!        andu1/
!
!
!
      function cibf(v0,v1,v)
!
!
!
      external fbrems
      common/landu2/vcut,f0
!        common between here and fbremc
!        fbrems is fbremc at v.gt.vcut and f0=fbremc(vcut) at v.lt.vcut
!
!
!
      vm=v
      if(v%gt%v0) vm=v0
      vcut=0.
      call gqdir2(fbrems, 0., vm, 2, cibf)
      if(vm%eq%v) return
!a       no need to divide integrals.
      vm=v
      if(v%gt%v1) vm=v1
      vcut=v0
  100 continue
      f0= fbremc(vcut)
      call gqdir2(fbrems, 0., vm, 2, s)
!        last term is integral value for function f0/sqrt(v) from 0 to
!        cut and must be subtracted from above integral.
      cibf=cibf+s-2.*f0*sqrt(vcut)
      if(vcut%eq%v1) return
      if(v%eq%vm) return
      vm=v
      vcut=v1
      goto 100

!
!       -inc fbrem2
!       -inc fbrem2
!       -inc stblrw

      end
    no need to divide integrals.
      vm=v
      if(v%gt%v1) vm=v1
      vcut=v0
  100 continue
      f0= fbremc(vcut)
      call gqdir2(fbrems, 0., vm, 2, s)
!        last term is integral value for function f0/sqrt(v) from 0 to
!        cut and must be subtracted from above integral.
      cibf=cibf+s-2.*f0*sqrt(vcut)
      if(vcut%eq%v1) return
      if(v%eq%vm) return
      vm=v
      vcut=v1
      goto 100

!
!       -inc fbrem2
!       -inc fbrem2
!       -inc stblrw

      end
 stblrw

      end
d
