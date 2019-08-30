
      subroutine ghstop(ipart,code,istop)
      integer *4 code,istop,ipart
c
c *** handling of stopping particles ***
c *** nve 18-may-1988 cern geneva ***
c
c called by : gheish
c origin : h.fesefeldt (routine calim 16-sep-1987)
c
c
c
c
c
c --- "ipart" changed to "kpart" in common /result/ due to clash ---
c --- with variable "ipart" in geant common ---
c
      common/result/xend,yend,zend,rca,rce,amas,nch,tof,px,py,pz,
     $              userw,intct,p,en,ek,amasq,deltn,itk,ntk,kpart,ind,
     $              lcalo,icel,sinl,cosl,sinp,cosp,
     $              xold,yold,zold,pold,pxold,pyold,pzold,
     $              xscat,yscat,zscat,pscat,pxscat,pyscat,pzscat
                    real nch,intct
c
c
c --- in case of energy deposition all the ekin will be deposited ---
      edep=ek
c
c --- update momentum vector and energies for stopping particle ---
      p=0.0
      px=0.0
      py=0.0
      pz=0.0
      en=abs(amas)
      ek=0.0
      getot=en
      gekin=ek
      istop=2
c
c
c *** select process for current particle ***
c
c --- skip exotic particles ---
      if (ipart .ge. 48) go to 9999
c
c --- look for particles with special treatment ---
      if (ipart .eq. 9) go to 90
      if (ipart .eq. 12) go to 120
      if (ipart .eq. 13) go to 130
      if (ipart .eq. 15) go to 150
      if (ipart .eq. 25) go to 250
c
c --- only deposit all kinetic energy for p and heavy fragments ---
      if (ipart .eq. 14) go to 140
      if (ipart .ge. 45) go to 140
c
c --- let all other particles decay ---
coff  call gdecay  !leave this to gismo
      code=5
      istop=1
      go to 9999
c
c --- pi- absorbed by nucleus ---
 90   continue
      call pimabs(nopt)
      code=16
      istop=1
      go to 9999
c
c --- k- absorbed by nucleus ---
 120  continue
      call kmabs(nopt)
      code=16
      istop=1
      go to 9999
c
c --- neutron captured by nucleus ---
 130  continue
      if (edep .ge. 1.e-9) go to 9999
      call captur(nopt)
      code=18
      istop=1
      go to 9999
c
c --- anti-proton ==> annihilation ---
 150  continue
      call pbanh(nopt)
      code=17
      istop=1
      go to 9999
c
c --- anti-neutron ==> annihilation ---
 250  continue
      call nbanh(nopt)
      code=17
      istop=1
      go to 9999
c
c --- p or heavy fragment ==> only deposit kinetic energy ---
 140  continue
      code=19
      istop=2
c
 9999 continue
c
      return
      end


