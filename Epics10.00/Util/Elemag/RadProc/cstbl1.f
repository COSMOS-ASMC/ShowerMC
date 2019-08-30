!     ****************************************************************
!     *                                                              *
!     * cstblbp1:  create sampling table for brems and pair          *
!     *            at low energy where no landau effect exists       *
!     *            all other correction included                     *
!     *                                                              *
!     *                                                              *
!     ********************** cstbl$1 *********************************
!
!       prepare:  fbrem ( -inc works for this end; fbrem1 module
!                              used)
!                 stblrw (prepare in #l.load etc to link it)
!                 make jpunch^=0, if want create data statement in
!                 disk //ft07f001
!                 table data is created in ftxxf001 if xx>0, which can
!                 be seen by stbldisp program.   jpunch and z and xx
!                 can be speicfied by tss conversation.
! *** see open below
!
!     lines between --  should be used in sampling routines
!     change dimension  bla(xx, yy) to bla(xx*yy) etc in that case.
! ------------------------- partical screening ----------------------
!           *** electron ***
!
!           ! egmin (gev)  for brems.
      data egmnl/.1e-3/
!
      parameter (eblmin=1.e-3, eblmax=50., nbtcl=40)
      dimension tcbl(nbtcl)
!
!             brem:   region a
!
      parameter (eela1=eblmin, eela2=eblmax,neela=20,
     * uela1=.05, uela2=1., nuela=20, duela=(uela2-uela1)/(nuela-1) )
      dimension bla(nuela, neela)
!
!             brem:  region b  (sqrt(u))
      parameter (eelb1=eblmin, eelb2=eblmax, neelb=neela,
     * uelb1=0.,  uelb2=.2236068, nuelb=20,
     * duelb=(uelb2-uelb1)/(nuelb-1) )
      dimension blb(nuelb, neelb)
!
!           ** gamma ***
!
!                 total cross-section
      parameter (eplmin=3.e-3, eplmax= 30., nptcl=16)
      dimension tcpl(nptcl)
!
!          pair: region a
      parameter ( egla1=eplmin,  egla2=eplmax, negla=16,
     *   ugla1=.05, ugla2=1., nugla=16,dugla=(ugla2-ugla1)/(nugla-1) )
      dimension pla(nugla, negla)
!
!            pair: region b     sqrt(log10(e)) sqrt(u)
      parameter (eglb1=eplmin, eglb2=eplmax, neglb=15,
     *  uglb1=0., uglb2=.2236068, nuglb=16,
     *  duglb=(uglb2-uglb1)/(nuglb-1))
      dimension plb(nuglb, neglb)
!
! -------------------------------------------------------------------
!
      dimension etcbl(nbtcl), etcbll(nbtcl)
      dimension ebla(neela), eblal(neela)
      dimension eblb(neelb), eblbl(neelb)
      dimension etcpl(nptcl), etcpll(nptcl)
      dimension epla(negla), eplal(negla)
      dimension eplb(neglb), eplbl(neglb)
!
      data eps/1.e-3/
!
      character matter*4, cap*100, media*20, mem*6, file*25
!
      common / bpcom /  al183z,e,ccz,emass,bcoef,fz ,z333
!
      dimension ua(100)
!
      external bremr  ,  pairr
!
      common/upsic/upsi,vmax
!
      jpunch = 0
      io=0
      write(7, *) 'enter  matter '
      read(*, *) matter
      if(matter .eq. 'pb') then
          z=82.
      elseif(matter .eq. 'fe') then
          z=26.
      elseif(matter .eq. 'w') then
          z=74.
      elseif(matter .eq. 'air') then
          z=7.37
      elseif(matter .eq. 'cu') then
          z=29.
      elseif(matter .eq. 'c') then
          z=6.
      elseif(matter .eq. 'h2o') then
          z=7.23
      else
          write(7, '('' enter z for matter='',a4)') matter
          read(*, *) z
      endif
!
      write(7,'('' matter='',a4, '' z='', f7.2)') matter, z
      write(media,'('' media='',a4,'' z='',f6.2)') matter, z
!
!
      dummy=zpart( z )
!
!
      debtcl=log10(eblmax/eblmin)/(nbtcl-1)
      deela=log10(eela2/eela1)/( neela-1)
      deelb=log10(eelb2/eelb1)/( neelb-1)
!
      deptcl=log10(eplmax/eplmin)/( nptcl-1)
      degla=log10(egla2/ egla1) /( negla-1)
      deglb=sqrt( log10(eglb2/ eglb1) )/( neglb-1)
!
          write(*, '(''      data debtcl/'',f9.6,''/'')') debtcl
          write(*, '(''      data deela/'',f9.6,''/'')') deela
          write(*, '(''      data deelb/'',f9.6,''/'')') deelb
!
!
          write(*, '(''      data deptcl/'',f9.6,''/'')') deptcl
          write(*, '(''      data degla/'',f9.6,''/'')') degla
          write(*, '(''      data deglb/'',f9.6,''/'')') deglb
!
      de1=10.**debtcl
      e=eblmin
      do   i=1, nbtcl
          etcbl(i)=e
          etcbll(i)=log10(e)
          vmin=egmnl/e
          vmax=vmaxv(e)
          call totcb(vmin,vmax,tcbl(i) )
          e=e*de1
      enddo
!
!
!
          el1=log10(eblmin)
          el2=log10(eblmax)
          write(*,'(''c        electron energy=''/(1hc,5x, 5g12.5))')
     *    etcbl
          write(*, '(''c       log10 step='', f8.4)') debtcl
          write(*, '(''c'')')
          write(*, '(''c      log10 of electron energy boundary'')')
          write(*, '(''      data  etcbl1/'',f9.6,''/'')') el1
          write(*, '(''c    *     ,etcbl2/'',f9.6,''/'')') el2
          write(*, '(''c'')')


       call kmkDataStm1(tcbl, nbtcl, 'tcbl', 'f7.4', 7)
!
!
!         region a
!
      de1=10.**deela
      e=eela1
       do   ie=1,neela
          ebla(ie)=e
          eblal(ie)=log10(e)
          vmin=egmnl/e
          vmax=vmaxv(e)
          call totcb(vmin,vmax,tcb)

          vl=vmin
          vr=vmax
          u=uela1
          do   iu=1,nuela -1
             upsi=u*tcb
             call kbchop(bremr,vl,vr,eps,v,j)
             if(j .le.0) then
                 write(0, *) ' cond j, e, u=', j, e, u
                 do vx = vmin, vmax, (vmax-vmin)/100.d0
                    call totcb(vx, vmax, tcbx)
                    write(0, *) vx, tcbx
                 enddo
                 write(0,*) 'u upsi, vl, vr', u, upsi,
     *           vl, vr, ' vmin=', vmin, ' tcb=', tcb
                 call totcb(vl,vmax,tcbleft)
                 call totcb(vr,vmax,tcbright)
                 write(0, *) ' tcb at vl, vr=', tcbleft,
     *           tcbright
            endif
             bla(iu,ie)=alog(v/vmin)/(1.-u)
             u=u+duela
          enddo
          bla(nuela,ie)=tcb/vmin/fbrem(vmin)
          e=e*de1
       enddo
!
      u=uela1
      do   i=1, nuela
         ua(i)=u
         u=u+duela
      enddo
!
!
      el1=log10(eela1)
      el2=log10(eela2)

          write(*, '(''c'')')
          write(*, '(''c       brem: region a'')')
          write(*,'(''c        electron energy=''/(1hc,5x, 5g12.5))')
     *    ebla
          write(*, '(''c       log10 step='', f8.4)') deela
          write(*, '(''c'')')
          write(*, '(''c      log10 of electron energy boundary'')')
          write(*, '(''      data  eela1l/'',f9.6,''/'')') el1
          write(*, '(''c    *     ,eela2l/'',f9.6,''/'')') el2
          write(*, '(''c'')')

!
!

          write(*, '(''c'')')
          write(*,'(''c            bla(iu,ie)=alog(v/vmin)/(1.-u)'')')
          write(*,'(''c            u='',f5.3,'' to '',f5.3,
     *    '' step='',f7.4)') uela1, uela2, duela
          write(*,'(''c            log10(e)='',f8.4, '' to '',f8.4,
     *    '' step='', f7.4)')  el1, el2, deela
          write(*, '(''c'')')


       call kmkDataStm1(bla, neela*nuela, 'bla', 'f8.4', 8)
!
!
!        region b
!
      e=eelb1
      de=10.**deelb
       do   ie=1,neelb
          eblb(ie)=e
          eblbl(ie)=log10(e)
          vmin=egmnl/e
          vmax=vmaxv(e)
          call totcb(vmin,vmax,tcb)
          vl=vmin
          vr=vmax-eps
          us=uelb1 + duelb
           do   iu=2,nuelb
              u=us**2
              upsi=u*tcb
              call kbchop(bremr,vl,vr,eps,v,j)
              if(j .le. 0) then
                  write(0,'('' e,u='',2g12.3)') e, u
              endif
              blb(iu,ie)=alog(v/vmin)
              us=us+duelb
           enddo
          blb(1, ie)=alog(vmax/vmin)
          e=e*de1
       enddo
!
      us=uelb1
       do   i=1, nuelb
         ua(i)=us
         us=us+duelb
       enddo
!
!
!
      el1=log10(eelb1)
      el2=log10(eelb2)

          write(*, '(''c'')')
          write(*, '(''c       brem: region b'')')
          write(*,'(''c        electron energy=''/(1hc,5x, 5g12.5))')
     *    eblb
          write(*, '(''c       log10 step='', f8.4)') deelb
          write(*, '(''c'')')
          write(*, '(''c      log10 of electron energy boundary'')')
          write(*, '(''      data  eelb1l/'',f9.6,''/'')') el1
          write(*, '(''c    *     ,eelb2l/'',f9.6,''/'')') el2
          write(*, '(''c'')')

!
!

          write(*, '(''c'')')
          write(*,'(''c            blb(iu,ie)=alog(v/vmin)/(1.-u)'')')
          write(*,'(''c            sqrt(u)='',f5.3,'' to '',f5.3,
     *    '' step='',f7.4)') uelb1, uelb2, duelb
          write(*,'(''c            log10(e)='',f8.4, '' to '',f8.4,
     *    '' step='', f7.4)')  el1, el2, deelb
          write(*, '(''c'')')


       call kmkDataStm1(blb, neelb*nuelb, 'blb', 'f8.4', 8)
!


!         pair
!
!
!        total cross section
!
      de1=10.**deptcl
      e=eplmin
      vmin=.5
       do   ie=1, nptcl
          etcpl(ie)=e
          etcpll(ie)=log10(e)
          vmax=vmaxv(e)
          call totcp(vmin,vmax,tcp)
          tcpl(ie)=tcp*2.
          e=e*de1
       enddo
!
!
      el1=log10(eplmin)
      el2=log10(eplmax)

          write(*, '(''c'')')
          write(*, '(''c       pair: total x-section'')')
          write(*,'(''c        gamm energy=''/(1hc,5x, 5g12.5))')
     *    etcpl
          write(*, '(''c       log10 step='', f8.4)') deptcl
          write(*, '(''c'')')
          write(*, '(''c      log10 of gamma energy boundary'')')
          write(*, '(''      data  etcpl1/'',f9.6,''/'')') el1
          write(*, '(''c    *     ,etcpl2/'',f9.6,''/'')') el2
          write(*, '(''c'')')


       call kmkDataStm1(tcpl, nptcl, 'tcpl', 'f8.4', 8)
!

!
!          region a
!
      vmin=.5
      de1=10.**degla
      e=egla1
       do   ie=1,negla
          epla(ie)=e
          eplal(ie)=log10(e)
          vmax=vmaxv(e)
          call totcp(vmin,vmax,tcp)
          vl=vmin
          vr=vmax
          u=ugla1
           do   iu=1,nugla-1
              upsi=u*tcp
              call kbchop(pairr,vl,vr,eps,v,j)
              if(j .le. 0) then
                  write(6,'('' e, u='',2g12.3)') e,u
              endif
              pla(iu,ie)=(v-.5)/(1.-u)
              u=u+dugla
           enddo
          pla(nugla,ie)=tcp/fpair(.5)
!
          e=e*de1
       enddo
!
      u=ugla1
       do   i=1, nugla
         ua(i)=u
         u=u+dugla
       enddo
!
!
!
!
!
      el1=alog10(egla1)
      el2=alog10(egla2)

           write(*, '(''c'')')
           write(*, '(''c       pair: region a'')')
           write(*, '(''c'')')
           write(*,'(''c         energy of gamma=''/(1hc,5x,5g12.3))')
     *     epla
           write(*,'(''c         log10 step='',f9.6)') degla
           write(*,'(''c'')')
           write(*,'(''c         log10 of eg boundary'')')
           write(*,'(''      data  egla1l/'',f9.6,''/'')') el1
           write(*,'(''c    *     ,egla2l/'',f9.6,''/'')') el2
           write(*, '(''c'')')

!

          write(*,'(''c'')')
          write(*,'(''c         pla(iu,ie)=(v-.5)/(1.-u)'')')
          write(*,'(''c         u='',f5.3,'' to'',f5.3,
     *    ''   step='',f7.4)')  ugla1, ugla2, dugla
          write(*,'(''c            log10(e)='',f8.4, '' to '',f8.4,
     *    '' step='', f7.4)')  el1, el2, degla
          write(*, '(''c'')')


       call kmkDataStm1(pla, nugla*negla, 'pla', 'f8.4', 8)
!

!     region b
!
      vmin=.5
      de1=deglb
      ex=0.
       do   ie=1,neglb
          e=10.**(ex**2) * eglb1
          eplb(ie)=e
          eplbl(ie)=log10(e)
          vmax=vmaxv(e)
          call totcp(vmin,vmax,tcp)
          vl=vmin
          vr=vmax
          us=uglb1 + duglb
           do   iu=2, nuglb
              u=us**2
              upsi=u*tcp
              call kbchop(pairr,vl,vr,eps,v,j)
              if(j .le. 0) then
                  write(6,'('' e, u='',2g12.3)') e,u
              endif
              plb(iu,ie)=v
              us=us+duglb
           enddo
          plb( 1,   ie)=vmax
!
          ex=ex+de1
       enddo
!
      us=uglb1
       do   i=1, nuglb
         ua(i)=us
         us=us+duglb
       enddo
!
!
!
!
!
      el1=alog10(eglb1)
      el2=alog10(eglb2)

           write(*, '(''c'')')
           write(*, '(''c'')')
           write(*, '(''c       pair: region b'')')
           write(*, '(''c'')')
           write(*,'(''c         energy of gamma=''/(1hc,5x,5g12.3))')
     *     eplb
           write(*,'(''c         sqrt(log10) step='',f9.6)') deglb
           write(*,'(''c'')')
           write(*,'(''c         use 0 for e-side at interpolation'')')
           write(*,'(''c         with sqrt(log10(e)-eglb1l)'')')
           write(*,'(''      data  eglb1l/'',f9.6,''/'')') el1
           write(*,'(''c'')')

!

          write(*,'(''c'')')
          write(*,'(''c         plb(iu,ie)=v'')')
          write(*,'(''c         sqrt(u)='',f5.3,'' to'',f5.3,
     *    ''   step='',f7.4)')  uglb1, uglb2, duglb


       call kmkDataStm1(plb, nuglb*neglb, 'plb', 'f8.4', 8)
!


      end
!
!
!
!     *****************
      subroutine totcb(vmin,vmax,ans)
!     *****************
!
!        integration of bremsung function from vmin to =max.
!       for integration ifovx and gquadt are employed.
!
!     external fbremv,fbrem,fpair
      external fpair
      external  fbrem
      ans=0.
      v2=vmax
    5 continue
      v1=v2/10.
      v1=max(v1,vmin)
      call gquadt(fbrem,v1,v2,ans1)
!       call k16pGaussLeg(fbrem, v1, v2, 16,  ans1)
      ans=ans+ans1
      if(v1 .eq. vmin) return
      v2=v1
      goto 5
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
!      call k16pGaussLeg(fpair, vmin, vt, 10,  ans1)
      call gquadt(fpair,vt,vmax,ans2)
!      call k16pGaussLeg(fpair, vt, vmax, 10,  ans2)
      ans=ans1+ans2
      return
      end
!
!
!     **************
      function bremr(v)
!     **************
!
!          used to solve total-cross-section * u = integral of
!          brem or pair function from min to v.
!
!
!
      common/upsic/upsi,vmax
      call totcb(v,vmax,ans)
! ////////////
!      write(*, *) ' v, ans=', v, ans, ' upsi=', upsi
!///////////////
   10 continue
      ans = max(0., ans)
      bremr=ans/upsi-1.
!       bremr = ans - upsi
      return
!
!     ************
      entry pairr(v)
!     ************
!
      call totcp(v,vmax,ans)
      goto 10
      end


r(v)
!     **************
!
!          used to solve total-cross-section * u = integral of
!          brem or pair function from min to v.
!
!
!
      common/upsic/upsi,vmax
      call totcb(v,vmax,ans)
! ////////////
!      write(*, *) ' v, ans=', v, ans, ' upsi=', upsi
!///////////////
   10 continue
      ans = max(0., ans)
      bremr=ans/upsi-1.
!       bremr = ans - upsi
      return
!
!     ************
      entry pairr(v)
!     ************
!
      call totcp(v,vmax,ans)
      goto 10
      end

