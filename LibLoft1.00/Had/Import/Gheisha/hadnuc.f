*cmz :          22/02/91  16.34.00  by  federico carminati
*cmz :  3.14/13 22/06/90  08.17.18  by  rene brun
*-- author :    rene brun   23/05/90
      subroutine hadnuc
     x(pvec,mass,atomic,ztomic,ipart,iok,stop,mprod,nprod,iprod,pprod)
c
      integer *4 stop
      real *4 pvec(3),ppvec,mass
      integer *4 mprod,nprod,iprod(mprod)
      real *4 pprod(3,mprod)
      integer *4 fail
c.
c.    ******************************************************************
c.    *                                                                *
c.    *           interface to the hadrin/nucrin package               *
c.    *                                                                *
c.    *                                                                *
c.    *       see paolo pedroni: infn/be-88/3   (6 april 88)           *
c.    *                                                                *
c.    *   "simulation of inelastic hadron collisions below 5 gev:"     *
c.    *          " a modification of the geant3 package"               *
c.    *                                                                *
c.    *                                                                *
c.    *       the pedroni's paper references the following articles:   *
c.    *                                                                *
c.    *     - k.hanssen and j.ranft, comp. phys. comm. 39, 37 (1986)   *
c.    *                                                                *
c.    *     - k.hanssen and j.ranft, comp. phys. comm. 39, 53 (1986)   *
c.    *                                                                *
c.    *                                                                *
c.    *                                                                *
c.    *    ==>called by : gheish                                       *
c.    *       author    r.brun , p.pedroni *********                   *
c.    *                                                                *
c.    ******************************************************************
c.
      common/gcflag/idebug,idemin,idemax,itest,idrun,idevt,ieorun
     +        ,ieotri,ievent,iswit(10),ifinit(20),nevent,nrndm(2)
c
      integer       idebug,idemin,idemax,itest,idrun,idevt,ieorun
     +        ,ieotri,ievent,iswit,ifinit,nevent,nrndm
c
      common/gcking/kcase,ngkine,gkin(5,100),tofd(100),iflgk(100)
      integer       kcase,ngkine ,iflgk
      real          gkin,tofd
c
      integer nmec,lmec,namec,nstep ,maxnst,ignext,inwvol,istop,maxmec
     + ,igauto,iekbin,ilosl, imull,ingoto,nldown,nlevin,nlvsav,istory
      real  vect,getot,gekin,vout,destep,destel,safety,sleng ,step
     + ,snext,sfield,tofg  ,gekrat,upwght
      parameter (maxmec=30)
      common/gctrak/vect(7),getot,gekin,vout(7),nmec,lmec(maxmec)
     + ,namec(maxmec),nstep ,maxnst,destep,destel,safety,sleng
     + ,step  ,snext ,sfield,tofg  ,gekrat,upwght,ignext,inwvol
     + ,istop ,igauto,iekbin, ilosl, imull,ingoto,nldown,nlevin
     + ,nlvsav,istory
c
      common/finuc/irn,itrn(60),cxrn(60),cyrn(60),czrn(60),
     +             elrn(60),plrn(60),tv
      common/finlsp/ir,itrr(20),cxr(20),cyr(20),czr(20),
     +              elr(20),plr(20)
c
c  ***** conversion tables between geant ipart and nucrin/hadrin it
c  --- igtonu(i)=nucrin  code corresponding to geant   code i ---
c  --- inutog(i)=geant   code corresponding to nucrin  code i ---
c
      dimension igtonu(48),inutog(23)
      data igtonu / 0, 0, 0, 0, 0, 0,23,13,14,25, 15,16, 8, 1, 2,24, 0,
     +17,21,22, 20, 0, 0, 0, 9,18,22*0/
     
      data inutog /14,15, 3, 2, 4, 4, 1,13,25, 5, 6,10, 8, 9,11,12,18,
     +26,16,21, 19,20, 7/
c.
c.    ------------------------------------------------------------------
c.
      if(ifinit(15).eq.0)then
         call acfsin
         call hadini
         call chanwx
         ifinit(15)=1
      endif
*
      iok=0
      ppvec=sqrt(pvec(1)**2+pvec(2)**2+pvec(3)**2)
      getot=sqrt(mass**2+ppvec**2)
      if (ppvec.gt.4.5) goto 999
      if(atomic.gt.1.5.and.atomic.lt.4.)goto 999
      if(ipart.gt.48)go to 999
      it=igtonu(ipart)
      if(it.eq.0) goto 999
      if(atomic.lt.1.5.and.it.le.22.and.it.ge.17)goto 999
      iok=1
      stop=1
      cx=pvec(1)/ppvec
      cy=pvec(2)/ppvec
      cz=pvec(3)/ppvec
      elab=getot
      if(atomic.lt.1.5) then
         plabo=vect(7)
         inuc=0
         call hadrinC(it,plabo,elab,cx,cy,cz,1,inuc)
         ngkine=ir
         do 10  i=1,ir
            fail=1301
            if(nprod.ge.mprod) go to 1313
            nprod=nprod+1
            pprod(1,nprod)=plr(i)*cxr(i)
            pprod(2,nprod)=plr(i)*cyr(i)
            pprod(3,nprod)=plr(i)*czr(i)
            if(itrr(i).lt.24) then
               iprod(nprod)=inutog(itrr(i))
            else
c -                             k0l/k0s choice from k0/k0bar
               call grndm(r,1)
               if(r.lt.0.5)then
                  iprod(nprod)=10.
               else
                  iprod(nprod)=16.
               endif
            endif
   10    continue
      else
         inuc=1
         call nucrinC(it,elab,cx,cy,cz,atomic,ztomic,0.)
**       destep=destep+tv
         ngkine=irn
         do 20  i=1,irn
            fail=1301
            if(nprod.ge.mprod) go to 1313
            nprod=nprod+1
            pprod(1,nprod)=plrn(i)*cxrn(i)
            pprod(2,nprod)=plrn(i)*cyrn(i)
            pprod(3,nprod)=plrn(i)*czrn(i)
            if(itrn(i).lt.24) then
               iprod(nprod)=inutog(itrn(i))
            else
c -                             k0l/k0s choice from k0/k0bar
               call grndm(r,1)
               if(r.lt.0.5)then
                  iprod(nprod)=10.
               else
                  iprod(nprod)=16.
               endif
            endif
   20    continue
      endif
*
  999 continue
      return
1313  continue
       write(6,*) 'fail=',fail,' in hadnuc'
       if(fail.eq.1301) then
        write(6,*) 'pprod is full; output will be truncated'
       endif !fail.eq.1301
      return
      end
