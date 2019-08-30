      function jamcomp(kf)

      implicit double precision(a-h, o-z)
      common/pydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/pydat4/chaf(500,2)
      character chaf*16
c     parameter(kfmx=100553)
      parameter(kfmx=35553)
      dimension kcode(kfmx)
      logical first
      save kcode
      save first
      data first/.true./

      if(first) then
        do i=1,kfmx
          kcode(i)=0
        end do
        first=.false.
        do kc=1,500
          kfa=kchg(kc,4)
          if(kfa.ge.1.and.kfa.le.kfmx) kcode(kfa)=kc
        end do
      endif

      if(kf.eq.0) then
        jamcomp=0
        return
      endif

      kfa=abs(kf)
      if(kfa.ge.1.and.kfa.le.kfmx) then
        if(mod(kfa/10,10).eq.0.and.kfa.lt.100000
     &     .and.mod(kfa/1000,10).gt.0) kfa=mod(kfa,10000)
        jamcomp=kcode(kfa)
      else
        jamcomp=jamcomp0(kf)
      endif


       end function
c***********************************************************************
 
      function jamcomp0(kf)
 
c...Compress the standard KF codes for use in mass and decay arrays;
c...also checks whether a given code actually is defined.

      implicit double precision(a-h, o-z)
      common/pydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/cpydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /pydat1/,/cpydat2/

c...Local arrays and saved data.
      dimension kford(100:500),kcord(101:500)
      save kford,kcord,nford,kflast,kclast
 
c...Whenever necessary reorder codes for faster search.
      if(mstu(20).eq.0) then
        nford=100
        kford(100)=0
        do i=101,500
          kford(i)=0
          kcord(i)=0
        end do

        do 120 i=101,500
          kfa=kchg(i,4)
          if(kfa.le.100) goto 120
          nford=nford+1
          do 100 i1=nford-1,0,-1
            if(kfa.ge.kford(i1)) goto 110
            kford(i1+1)=kford(i1)
            kcord(i1+1)=kcord(i1)
  100     continue
  110     kford(i1+1)=kfa
          kcord(i1+1)=i
  120   continue
        mstu(20)=1
        kflast=0
        kclast=0
      endif
 
c...Fast action if same code as in latest call.
      if(kf.eq.kflast) then
        jamcomp0=kclast
        return
      endif
 
c...Starting values. Remove internal diquark flags.
      jamcomp0=0
      kfa=iabs(kf)
      if(mod(kfa/10,10).eq.0.and.kfa.lt.100000
     &     .and.mod(kfa/1000,10).gt.0) kfa=mod(kfa,10000)
 
c...Simple cases: direct translation.
      if(kfa.gt.kford(nford)) then
      elseif(kfa.le.100) then
        jamcomp0=kfa
 
c...Else binary search.
      else
        imin=100
        imax=nford+1
  130   iavg=(imin+imax)/2
        if(kford(iavg).gt.kfa) then
          imax=iavg
          if(imax.gt.imin+1) goto 130
        elseif(kford(iavg).lt.kfa) then
          imin=iavg
          if(imax.gt.imin+1) goto 130
        else
          jamcomp0=kcord(iavg)
        endif
      endif
 
c...Check if antiparticle allowed.
      if(jamcomp0.ne.0.and.kf.lt.0) then
        if(kchg(jamcomp0,3).eq.0) jamcomp0=0
      endif
 
c...Save codes for possible future fast action.
      kflast=kf
      kclast=jamcomp0
 
       end function
