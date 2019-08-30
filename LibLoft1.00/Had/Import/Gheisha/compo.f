*cmz :  3.14/16 28/09/90  09.25.40  by  nick van eijndhoven (cern)
*-- author :
      subroutine compo
      integer *4 maxmix
      parameter (maxmix=5)
c
c *** get parameters for the target nucleus ***
c *** nve 11-mar-1988 cern geneva ***
c
c origin : h.fesefeldt (07-dec-1984)
c
      common/curpar/weight(10),ddeltn,ifile,irun,nevt,nevent,shflag,
     *              ithst,ittot,itlst,ifrnd,tofcut,cmom(5),ceng(5),
     *              rs,s,enp(10),np,nm,nn,nr,no,nz,ipa(200),
     *              atno2,zno2
c
c
      integer *4 nmix,imix(maxmix)
      real *4 a,z,wmix(maxmix),zmix(maxmix),amix(maxmix)
      character *20 name,namep,c20
      save name,a,z,nmix,imix,wmix
      save amix,zmix
      integer *4 nmixp,imixp,imatp
      real *4 ap,zp
      integer *4 ndum
      real *4 rndm(1)
      real *4 dum(10)
      integer *4 comp
      integer *4 fail
      data a/0/
      data z/0/
      data nmix/0/
c
      fail=1301
      if(a.eq.0) go to 1313
c
c
c --- check for compound ---
      kk=nmix
      if (kk .ge. 2) go to 10
c
c --- elements ---
      atno2=a
      zno2=z
      go to 9999
c
c --- compounds ===> select nucleus ---
 10   continue
c
      sum=0.0
      do 11 i=1,kk
      ai=amix(i)
      zi=zmix(i)
      wi=wmix(i)
      sum=sum+wi/ai
 11   continue
      call grndm(rndm,1)
      test1=rndm(1)*sum
c
      test2=0.0
      do 12 i=1,kk
      jcompo=i
      ai=amix(i)
      zi=zmix(i)
      wi=wmix(i)
      test2=test2+wi/ai
      if (test2 .gt. test1) go to 20
 12   continue
c
 20   continue
      atno2=ai
      zno2=zi
c
 9999 continue
      return
1313  continue
       write(6,*) 'fail=',fail,' in compo'
       if(fail.eq.1301) then
        write(6,*) 'compo was not inilized'
       endif !fail.eq.1301
      return
c
      entry compoi(imatp,namep,ap,zp,nmixp,imixp,wmixp)
       imat=imatp
       name=namep
       a=ap
       z=zp
       nmix=nmixp
       call ucopy(imixp,imix,nmixp)
       call ucopy(wmixp,wmix,nmixp)
       do 1000 i=1,nmix
        comp=imix(i)
        ndum=10
        call gfmate(comp,c20,amix(i),zmix(i),dum,dum,dum,ndum,dum,dum)
1000   continue
      return
      end
