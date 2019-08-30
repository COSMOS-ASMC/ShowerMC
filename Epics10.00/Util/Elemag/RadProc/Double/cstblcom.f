!     ****************************************************************
!     *                                                              *
!     * cstblcom: create sampling table for compton scattering       *
!     *                                                              *
!     ****************************************************************
!
!     prepare:  fcomp( -inc works for this end)
!               fix material by call zpart(z,a,rho)
!               give jpunch^=0 if want to produce data stmt in disk
!        to use gd, use gdstbl jcl in #tape.cntl
!
!      material   z     a     rho     t0 (g/cm2)  t0 (cm)    ec
!         pb      82   207.2  11.35    5.82        0.513     6.7
!         cu      29    63.5   8.94   12.7        1.42     16+x
!         fe      26    55.85  7.86   13.7         1.74     18+x
!          w      74   183.92 19.3     6.275       .325     8.08
!
!
!      lines between ---  should be included in sampling routien
!
!     -----------------------------------------------------------
!                              negc=log10(egc2/egc1)/degcl + 1(or 2)
      parameter (degcl=1./4., egc1=.1e-3, egc2=100.e-3, negc=13,
     * dugc=.05,
     * ugc1=0., ugc2=1., nugc=(ugc2-ugc1+.0001)/dugc+1)
      dimension tcc(negc)
!
!     -----------------------------------------------------------
!
      data stepf/3./,eps/1.e-3/
!
      dimension uec1(nugc, negc), uec2(nugc,negc), egm(negc)
!
*     dimension ua(nugc)
!
      common/ccomp/e,upsi
!
      common /$compt/  emass, cconst, x0g, x0cm
!
      external compr
!
      jpunch=1
!
      call zpart(26., 55.85,   7.86)
!
      lmx=6
      if(jpunch .ne. 0) lmx=7
!
      de=10.**degcl
      e=egc1
       do   i=1, negc
          egm(i)=e
          vmin=vminc(e, emass)
          ans=ficomp(e, emass, vmin)
          tcc(i)=ans
          u=dugc
          v=1.-0.1*(1.-vmin)
          vmax=1.
!
           do   j=2,nugc-1
              upsi=u*ans
              call rutbnr(compr,v,eps,vmin,vmax, stepf,v,k)
              uec1(j,i)= alog(v)
              uec2(j,i)= v
              vmax=v
              v=(vmin-v)/(1.-u)*.1+v
              if(v%le%vmin) v=v+eps
              u=u+dugc
           enddo
          uec1(nugc,i)=alog(vmin)
          uec2(nugc,i)=vmin
          uec1(1,i)=0.
          uec2(1,i)=1.
          e=e*de
       enddo
!            gd  write
*     u=0.
*     do 4900 i=1, nugc
*        ua(i)=u
*        u=u+dugc
*4900 continue
!
*     call stblgd(' compton:  alog(v); v=eg''/eg!',
*    * ua, uec1, nugc, negc, '  uec1    ',
*    * 12.5, 2. )
!
*     call stblgd(' comprton: v=eg''/eg!',
*    * ua, uec2, nugc, negc, ' uec2    ',
*    * 12.5, 10.)
!
!
!
!
!
       do   l=6, lmx
           egc1l=alog10(egc1)
           egc2l=alog10(egc2)
           write(l,'(''c'')')
           write(l,'(''c         energy of gamma=''
     *     /(1hc, 5x, 5g12.4))') egm
           write(l,'(''c         log10 step='',f9.6)') degcl
           write(l,'(''c         log10 gamma energy boundary'')')
           write(l,'(''      data egc1l/'',f9.6,''/, egc2l/'',
     *     f9.6,''/'')') egc1l, egc2l
           write(l,'(''c'')')
       enddo

       do   l=6, lmx
            write(l,'(''c'')')
            write(l,'(''c           uec1(j,i)=alog(v)'')')
            write(l,'(''c           from u='',f7.4, '' to'',f7.4,
     *      '' step='',f7.4)') ugc1,ugc2, dugc
            write(l,'(''c           from log10(e)='',f9.6,'' to'',f9.6,
     *      '' step='',f8.4)') egc1l, egc2l, degcl
            write(l,'(''c           dim. of u is'',i3, ''  e is'',i3)')
     *      nugc, negc
            write(l,'(''c'')')
            write(l,'(''      dimension uec1('',i4,'')'')') nugc*negc
            write(l,'(''c'')')
       enddo
      call mkdt(' uec1 ',uec1, 1,  nugc*negc,0., 0,  jpunch)

       do   l=6, lmx
            write(l,'(''c'')')
            write(l,'(''c           uec2(j,i)=v'')')
            write(l,'(''c           from u='',f7.4, '' to'',f7.4,
     *      '' step='',f7.4)') ugc1,ugc2, dugc
            write(l,'(''c           from log10(e)='',f9.6,'' to '',f9.6
     *      '' step='',f8.4)') egc1l, egc2l, degcl
            write(l,'(''            dim. of u and e='',2i4)')
     *      nugc, negc
            write(l,'(''c'')')
            write(l,'(''      dimension uec2('',i4,'')'')') nugc*negc
            write(l,'(''c'')')
       enddo
      call mkdt(' uec2 ',uec2, 1,  nugc*negc,'f7.5,  ', 7,jpunch)
!
!
!
!
*     call ggdtm
      stop
      end
      function compr(v)
      common /$compt/  emass, cconst, x0g, x0cm
      common/ccomp/e,upsi
!
      compr=ficomp(e, emass, v)/upsi-1.
      return
      end
!
!       -inc fcomp
       -inc fcomp
st, x0g, x0cm
      common/ccomp/e,upsi
!
      compr=ficomp(e, emass, v)/upsi-1.
      return
      end
!
!       -inc fcomp
       -inc fcomp
