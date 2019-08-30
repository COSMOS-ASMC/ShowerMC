*cmz :  3.14/16 29/06/89  14.19.16  by  nick van eijndhoven (cern)
*-- author :
      real function atomas(a,z)
c
c *** determination of the atomic mass ***
c *** nve 19-may-1988 cern geneva ***
c
c origin : h.fesefeldt (02-dec-1986)
c
      common/consts/ pi,twpi,pibtw,mp,mpi,mmu,mel,mkch,mk0,smp,smpi,
     $               smu,ct,ctkch,ctk0,
     $               ml0,msp,ms0,msm,mx0,mxm,ctl0,ctsp,ctsm,ctx0,ctxm,
     $               rmass(35),rcharg(35)
c
                     real mp,mpi,mmu,mel,mkch,mk0,
     *                    ml0,msp,ms0,msm,mx0,mxm
c
c
      double precision aa,zz,mass
c
c --- get atomic (= electrons incl.) masses (in mev) from rmass array ---
c --- electron ---
      rmel=rmass(4)*1000.
c --- proton ---
      rmp=rmass(14)*1000.
c --- neutron ---
      rmn=rmass(16)*1000.
c --- deuteron ---
      rmd=rmass(30)*1000.+rmel
c --- alpha ---
      rma=rmass(32)*1000.+2.*rmel
c
      atomas = 0.
      aa = a * 1.d0
      zz = z * 1.d0
      ia = ifix(a + 0.5)
      if(ia.lt.1) return
      iz = ifix(z + 0.5)
      if(iz.lt.0) return
      if(iz.gt.ia) return
      if(ia.gt.4) goto 50
      mass=0.d0
      goto (10,20,50,40),ia
   10 if(iz.eq.0) mass=rmn
      if(iz.eq.1) mass=rmp+rmel
      goto 60
   20 if(iz.ne.1) goto 50
      mass=rmd
      goto 60
   40 if(iz.ne.2) goto 50
      mass=rma
      goto 60
   50 mass=(aa-zz)*rmn + zz*rmp +zz*rmel - 15.67*aa
     *     + 17.23*(aa**0.6666667) + 93.15*((aa/2.-zz)**2)/aa
     *     +0.6984523*zz**2/(aa**0.3333333)
      ipp=mod(ia-iz,2)
      izz=mod(iz,2)
      if(ipp.ne.izz) goto 60
      mass = mass + (ipp+izz- 1)*12.00*(aa**(-0.5))
   60 atomas = mass*0.001
      return
      end
