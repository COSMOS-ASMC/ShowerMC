      subroutine NumberModel(cmodel,model)
      character cmodel*21
      n=index(cmodel,' ')-1
      if(cmodel(1:n).eq.'qgsjet' )   model=2
      if(cmodel(1:n).eq.'gheisha')   model=3
c      if(cmodel(1:n).eq.'gheisha' )stop'wrong choice!!!!'  !    model=4
      if(cmodel(1:n).eq.'pythia' )stop'wrong choice!!!!'  !    model=4
      if(cmodel(1:n).eq.'hijing' )stop'wrong choice!!!!'  !    model=5
      if(cmodel(1:n).eq.'sibyll' )   model=6
      if(cmodel(1:n).eq.'IIqgsjet' ) model=7
      end

      subroutine IniModel(model)
      if(model.eq.2)call IniQGSjet
      if(model.eq.3)call IniGheisha
cxxx      if(model.eq.4)call IniPythia
cxxx      if(model.eq.5)call IniHijing
      if(model.eq.6)call IniSibyll
      if(model.eq.7)call IniQGSJetII
      end

      subroutine IniEvtModel
      include 'epos.inc'
      if(model.eq.2)call IniEvtQGS
      if(model.eq.3)call IniEvtGhe
cxxx      if(model.eq.4)then
cxxx        engysave=engy
cxxx        if(engy.lt.egymin)engy=egymin
cxxx        call IniEvtPyt
cxxx        engy=engysave
cxxx      endif
c      if(model.eq.5)then
c        engysave=engy
c        if(engy.lt.egymin)engy=egymin
c        call IniEvtHij
c        engy=engysave
c      endif
      if(model.eq.6)call IniEvtSib
      if(model.eq.7)call IniEvtQGSII
      end

      subroutine emsaaaModel(model,id,iret)
        if(model.eq.2)then
          if(id.eq.0)call IniEvtQGS
          call emsqgs(iret)
        endif
        if(model.eq.3)call emsghe(iret)
cxxx    if(model.eq.4)call emspyt(iret)
c       if(model.eq.5)call emshij(iret)
        if(model.eq.6)call emssib(iret)
        if(model.eq.7)then
          if(id.eq.0)call IniEvtQGSII
          call emsqgsII(iret)
        endif
      end

      function crseModel(model,ekin,maproj,matarg,idtarg)
        if(model.eq.2)crseModel=qgscrse(ekin,maproj,matarg,idtarg)
cxxx        if(model.eq.4)??????????????????
c        if(model.eq.5)crseModel=hijcrse(ekin,maproj,matarg,idtarg)
        if(model.eq.6)crseModel=sibcrse(ekin,maproj,matarg,idtarg)
        if(model.eq.7)crseModel=qgsIIcrse(ekin,maproj,matarg,idtarg)
      end

      subroutine crseaaModel(sigt,sigi,sigc,sige)
      include 'epos.inc'
      if(idtarg.eq.0)then
        sigt=0.
        sigc=0.
        sigi=0.
        sige=0.
        do k=1,3
          matarg=int(airanxs(k))
          call crseaaModel0(xsigt,xsigi,xsigc,xsige)
          sigt=sigt+airwnxs(k)*xsigt
          sigi=sigi+airwnxs(k)*xsigi
          sigc=sigc+airwnxs(k)*xsigc
          sige=sige+airwnxs(k)*xsige
        enddo
      else
        call crseaaModel0(sigt,sigi,sigc,sige)
      endif
      end

      subroutine crseaaModel0(sigt,sigi,sigc,sige)
      include 'epos.inc'

      double precision GTOT,GPROD,GABS,GDD,GQEL,GCOH
      if(model.eq.2)then
        NITER=5000
        if(idtarg.eq.0)then
          e0=dble(elab)
          icp=idtrafo('nxs','qgs',idproj)
          call xxaini(e0,icp,maproj,matarg)
        endif
        CALL CROSSC(NITER,GTOT,GPROD,GABS,GDD,GQEL,GCOH)
        sigt=sngl(GTOT)
        sigi=sngl(GPROD)
        sigc=sngl(GABS)
        sige=sigt-sigi
      elseif(model.eq.3)then
        idtar=idtarg
        if(idtarg.eq.0)idtar=1120
        call ghecrse(ekin,idproj,idtar,latarg,matarg,sigi,sige)
        sigt=sigi+sigel
        sigc=sigi
      elseif(model.eq.6.and.maproj.eq.1)then
        K=1
        if(iclpro.eq.1)then
          K=2
        elseif(iclpro.eq.3)then
          K=3
        endif
        CALL SIB_SIGMA_HP(K,engy,SSIG,dum0,dum1,dum2,SLOPE,RHO)
        CALL GLAUBER(matarg,SSIG,SLOPE,RHO,sigt,sige,sigqel)
        sigi=sigt-sigqel
        sigc=sigi
      elseif(model.eq.7)then
        NITER=5000
        if(idtarg.eq.0)then
          e0=dble(elab)
          icp=idtrafo('nxs','qgs',idproj)
          call qgini(e0,icp,maproj,matarg)
        endif
        CALL qgcrossc(NITER,GTOT,GPROD,GABS,GDD,GQEL,GCOH)
        sigt=sngl(GTOT)
        sigi=sngl(GPROD)
        sigc=sngl(GABS)
        sige=sigt-sigi
      else
        sigt=0.
        sigi=0.
        sigel=0.
        sigc=0.
      endif
      end


      subroutine m2XXFZ( a,b)
      double precision a,b
        CALL XXFZ(a,b)
      end

      subroutine m3SIGMA(ek,idpro,idtar,latar,matar,sigi,sige)
      call ghecrse(ek,idpro,idtar,latar,matar,sigi,sige)
      end

      subroutine m6SIGMA(icl,engy,stot,sela,sine,sdifr,slela,Rho)
        if(icl.eq.1)then
          L=2
        elseif(icl.eq.2)then
          L=1
        else
          L=3
        endif
        call SIB_SIGMA_HP(L,engy,stot,sela,sine,sdifr,slela,Rho)
      end


      subroutine m7SIGMA(stot,scut,sine,slela)
      double precision GzZ0(5),pi,bm,amws
      common /qgarr1/  ia(2),icz,icp
      common /qgarr6/  pi,bm,amws
      ia2save=ia(2)
      bmsave=bm
      ia(2)=1
      call qgfz(0.d0,gzz0,0,0)
      scut=sngl(gzz0(2))/2.   !cut pomerons cross-section
      stot=sngl(gzz0(1))                     !tot cross-section
      sine=sngl(gzz0(2)+gzz0(3)+gzz0(4))/2. !inelastic cross section
      slela=sngl(gzz0(5))
      ia(2)=ia2save
      bm=bmsave
      end
