c***********************************************************************
c                                                                      *
c        PART  : BUU part                                              *
c                                                                      *
c   List of subprograms in rough order of relevance with main purpose  *
c      (S = subroutine, F = function, B = block data, E = entry)       *
c                                                                      *
c  s  jambuuin to initialize BUU part                                  *
c  s  jambuud  to calculate nuclear density                            *
c  e  jambuudi to initialize gaussian for calculation of density       *
c  s  jammxdns to calculate maximum density                            *
c  s  jambuue  to calculate total energy per nucl.                     *
c  e  jambuur  to calculate rms.                                       *
c  s  jambuuf  to determine nuclear force acting on the hadrons        *
c  s  jambuuci to  calculate the initial value of coulomb potential    *
c  s  jambuucb to calculate the coulomb potential on the boundary      *
c  s  jambuucl to calculate the coulomb potential                      *
c                                                                      *
c                                                                      *
c***********************************************************************
c***********************************************************************

c***********************************************************************

      subroutine jambuuin(msel)

c...Purpose: to initialize BUU part.

      implicit double precision(a-h, o-z)
      include 'jam2.inc'

c...msel=0: initialize over all event.
c...msel=1: initialize for each event.

      if(msel.eq.0) then

c...Initialize density part.
        call jambuudi
      endif

      if(mstc(6).ge.2) then

c...Soft e.o.s.
        if(mod(mstc(6),10).eq.2) then
          pard(101)=-0.3581d0
          pard(102)=0.3048d0
          pard(103)=7d0/6d0
        else if(mod(mstc(6),10).eq.4) then
c...Hard e.o.s.
          pard(101)=-0.124d0
          pard(102)=0.0705d0
          pard(103)=2d0
        else
          call jamerrm(30,0,'(jambuuin:)invalid mstc6')
        endif

c...Initialize Coulomb part.
      if(mstc(101).eq.1) call jambuuci(msel)

      endif

      end

c***********************************************************************

      subroutine jambuud

c...Purpose: to calculate nuclear density from spatial distribution of
c...testparticles and average momentum in spacial cell 125*125 points
c...gaussian smearing width of gaussian is 1.0 fm, cutoff=sqrt(5.0) fm.

      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
c...Gauss smearing.
      common/pq/cm(130,130),dm(130,350),minner,mouter,moutep

      ip=minner
      iq=mouter

      do iz = -maxz,maxz
      do iy = -maxx,maxx
      do ix = -maxx,maxx
        rhob(ix,iy,iz)=0.0d0
        rhoz(ix,iy,iz)=0.0d0
        avp(1,ix,iy,iz)=0.0d0
        avp(2,ix,iy,iz)=0.0d0
        avp(3,ix,iy,iz)=0.0d0
        avp(4,ix,iy,iz)=1.0d0
      end do
      end do
      end do

c...Loop over all particles.
      do 100 i=1,nv

      if(k(1,i).gt.10) goto 100
      if(pard(1)-r(4,i).lt.0.0d0) goto 100
c     if(pard(1)-r(5,i).lt.0.0d0) goto 100

      ix=nint(r(1,i))
      iy=nint(r(2,i))
      iz=nint(r(3,i))
      if(abs(ix).gt.maxx.or.abs(iy).gt.maxx.or.abs(iz).gt.maxz) goto 100

      kx=nint(dble(2*ip+1)*(r(1,i)-dble(ix)))
      if(abs(kx) .eq. ip+1) kx = kx/abs(kx) * ip

      ky=nint(dble(2*ip+1)*(r(2,i)-dble(iy)))
      if(abs(ky) .eq. ip+1) ky = ky/abs(ky) * ip

      kz=nint(dble(2*ip+1)*(r(3,i)-dble(iz)))
      if(abs(kz) .eq. ip+1) kz = kz/abs(kz) * ip

      ic=1+(kz+ip)+(ky+ip)*(2*ip+1)+(kx+ip)*(2*ip+1)**2
      ib=0

c...Loop for Gauss smearing.
      do 110 jx=ix-iq,ix+iq
      do 110 jy=iy-iq,iy+iq
      do 110 jz=iz-iq,iz+iq
        ib=ib+1
        if(cm(ic,ib).gt.0.0.and.
     &    abs(jx).le.maxx.and.abs(jy).le.maxx.and.abs(jz).le.maxz) then
 
c...baryon density.
c       kf=k(2,i)
c       id=kchg(jamcomp(kf),5)

c....Charge denstiy.
        kch=jamchge(k(2,i))
        if(kch.ne.0) then
          rhoz(jx,jy,jz)=rhoz(jx,jy,jz)+cm(ic,ib)*kch/3
        endif

        ibar=k(9,i)/3
        if(abs(ibar).eq.1) then
          rhob(jx,jy,jz)=rhob(jx,jy,jz)+cm(ic,ib)
          avp(1,jx,jy,jz)=avp(1,jx,jy,jz)+p(1,i)*cm(ic,ib)
          avp(2,jx,jy,jz)=avp(2,jx,jy,jz)+p(2,i)*cm(ic,ib)
          avp(3,jx,jy,jz)=avp(3,jx,jy,jz)+p(3,i)*cm(ic,ib)
          avp(4,jx,jy,jz)=avp(4,jx,jy,jz)+p(4,i)*cm(ic,ib)
        endif

        endif
  110 continue
  100 continue  ! end loop over all particles.
 
 
c...Calculate momentum and gamma factor of each cell.
      do 200 iz = -maxz,maxz
      do 210 iy = -maxx,maxx
      do 220 ix = -maxx,maxx

c...Density of all baryons.
        s=avp(4,ix,iy,iz)**2
     $-avp(1,ix,iy,iz)**2-avp(2,ix,iy,iz)**2-avp(3,ix,iy,iz)**2
        if(s.le.0.0d0) goto 220
        gam1=avp(4,ix,iy,iz)/sqrt(s)
        rhob(ix,iy,iz)=rhob(ix,iy,iz)/gam1

  220 continue
  210 continue
  200 continue

c...Calculate Coulomb potential.
      if(mstc(101).ge.1) then
        if((mstd(23)/mstc(102))*mstc(102).eq.mstd(23)) call jambuucb
        call jambuucl
      endif
 
      return
 
c***********************************************************************
                                                                        
      entry jambuudi

c....Initialize gaussian for calculation of density.
c     ip ; inner mesh for dens
c     iq ; outer mesh for dens
c     ir ; outer mesh for pauli and average momentum

      minner=2
      mouter=2
      moutep=3

      ip=2
      iq=2
      ir=3
 
      gauss=1.0d0
      dcut=5.0d0
 
      rd=1.0d0/dble(ip*2+1)
      do 101 ix=-ip,ip
      do 101 iy=-ip,ip
      do 101 iz=-ip,ip
      ic=1+(iz+ip)+(iy+ip)*(2*ip+1)+(ix+ip)*(2*ip+1)**2
      xi=dble(ix)*rd
      yi=dble(iy)*rd
      zi=dble(iz)*rd
      ib=0
      sek=0.0d0
      do 201 jx=-iq,iq
      do 201 jy=-iq,iq
      do 201 jz=-iq,iq
      ib=ib+1
      xj=dble(jx)
      yj=dble(jy)
      zj=dble(jz)
      rsqr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      if(rsqr.gt.dcut) then
      cm(ic,ib)=0.0d0
      else
      cm(ic,ib)=exp(-rsqr/2.d0/gauss**2)
      end if
      sek=sek+cm(ic,ib)
  201 continue
      ib=0
      do 301 jx=-iq,iq
      do 301 jy=-iq,iq
      do 301 jz=-iq,iq
      ib=ib+1
      cm(ic,ib)=cm(ic,ib)/sek/dble(mstc(5))
  301 continue
  101 continue
 
      do 121 ix=-ip,ip
      do 121 iy=-ip,ip
      do 121 iz=-ip,ip
      ic=1+(iz+ip)+(iy+ip)*(2*ip+1)+(ix+ip)*(2*ip+1)**2
      ie=0
      do 221 jx=-ir,ir
      do 221 jy=-ir,ir
      do 221 jz=-ir,ir
      ie=ie+1
      sek=0.0d0
      do 231 kx=jx-1,jx+1
      if(abs(kx).gt.iq) goto 231
      do 232 ky=jy-1,jy+1
      if(abs(ky).gt.iq) goto 232
      do 233 kz=jz-1,jz+1
      if(abs(kz).gt.iq) goto 233
      ib=1+(kz+iq)+(ky+iq)*(2*iq+1)+(kx+iq)*(2*iq+1)**2
      sek=sek+cm(ic,ib)
  233 continue
  232 continue
  231 continue
      dm(ic,ie)=sek
c     write(6,*) ic,ie,sek
  221 continue
  121 continue

      end

c***********************************************************************

      subroutine jammxdns(rhomx,rhoave,nvol)

c....Calculate maximum density.
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'

      rr1=pard(50)
      rr2=pard(40)
      rhomx=0.0d0
      rhoave=0.0d0
      nvol=0
      do iz = -maxz,maxz
      do iy = -maxx,maxx
      do ix = -maxx,maxx
       dr=sqrt(dble(ix)**2+dble(iy)**2+dble(iz)**2)
       if(dr.gt.rr1) cycle
       if(dr.gt.rr2) cycle
       rhotmp=rhob(ix,iy,iz)
       rhoave=rhoave+rhotmp
       rhomx=max(rhomx,rhotmp)
       nvol=nvol+1
      enddo
      enddo
      enddo

      if(nvol.gt.0) rhoave=rhoave/nvol

      end

c***********************************************************************

      subroutine jambuue(ekin,epot,etot)

c...Calculate total energy/A  with a fixed Eulerian method.
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'

c...Coff. of e.o.s.
      alp=pard(101)
      bet=pard(102)
      gam=pard(103)
      rho0=parc(21)

      ekin=0.0d0
      epot=0.0d0

c....Loop over all ensembles and particles.
      do 100 i=1,nv
        if(k(1,i).gt.10) goto 100
        ekin=ekin+sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)

        ix=nint(r(1,i))
        if(abs(ix).gt.maxx) goto 100
        iy=nint(r(2,i))
        if(abs(iy).gt.maxx) goto 100
        iz=nint(r(3,i))
        if(abs(iz).gt.maxz) goto 100

        if(mstc(101).ge.1) epot=epot+jamchge(k(2,i))/3*cop(ix,iy,iz)
        if(abs(k(9,i)).ne.3) goto 100
c       id=kchg(jamcomp(k(2,i),5)

        dens=rhob(ix,iy,iz)
cc      epot=epot+facv(1,inuc)*(alp/rho0*dens+bet*(dens/rho0)**gam)
        epot=epot+alp/rho0*dens+bet*(dens/rho0)**gam
 100  continue

      ekin=ekin/(mstc(5)*mstd(11))
      epot=epot/(mstc(5)*mstd(11))
      etot=ekin+epot

      return

c***********************************************************************

      entry jambuur(rms)

c...Calculate rms.

      rms=0.0d0
      do iz = -maxz,maxz
      do iy = -maxx,maxx
      do ix = -maxx,maxx
        rrt=rhob(ix,iy,iz)
        rms=rms+rrt*dble(ix**2+iy**2+iz**2)
      end do
      end do
      end do

      rms=sqrt(rms/mstc(5))

      end

c***********************************************************************

      subroutine jambuuf

c....Determine nuclear force acting on the hadrons.
      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'

c     parameter         (udel0=-0.03)
      parameter(dx=2.0d0, dy=2.0d0, dz=2.0d0)

c...Cof. of e.o.s.
      alp=pard(101)
      bet=pard(102)
      gam=pard(103)
      rho0=parc(21)

c...Cof. of skyrm potential.
      cof1=alp/rho0
      cof2=bet/rho0**gam

      do 100 i=1,nv
        if(k(1,i).gt.10) goto 100

        force(1,i)=0.0d0
        force(2,i)=0.0d0
        force(3,i)=0.0d0

        if(pard(1)-r(4,i).lt.0.0d0) goto 100
c       if(pard(1)-r(5,i).lt.0.0d0) goto 100

        ix=nint(r(1,i))
        iy=nint(r(2,i))
        iz=nint(r(3,i))
        if(abs(ix).ge.maxx.or.abs(iy).ge.maxx.or.abs(iz).ge.maxz)
     $  goto 100

c...Coulomb force.
        ich=jamchge(k(2,i))/3
        if(mstc(102).ge.1.and.ich.ne.0) then
        force(1,i)=-ich*(cop(ix+1,iy,iz)-cop(ix-1,iy,iz))/dx
        force(2,i)=-ich*(cop(ix,iy+1,iz)-cop(ix,iy-1,iz))/dy
        force(3,i)=-ich*(cop(ix,iy,iz+1)-cop(ix,iy,iz-1))/dz
        endif

c...Only baryons feel nuclear potential.
        if(abs(k(9,i)).ne.3) goto 100

        tdenxp = rhob(ix+1,iy,iz)
        tdenxm = rhob(ix-1,iy,iz)
        tdenyp = rhob(ix,iy+1,iz)
        tdenym = rhob(ix,iy-1,iz)
        tdenzp = rhob(ix,iy,iz+1)
        tdenzm = rhob(ix,iy,iz-1)

        rxderv = (tdenxp-tdenxm)/dx
        ryderv = (tdenyp-tdenym)/dy
        rzderv = (tdenzp-tdenzm)/dz

        dxderv = (tdenxp**gam-tdenxm**gam)/dx
        dyderv = (tdenyp**gam-tdenym**gam)/dy
        dzderv = (tdenzp**gam-tdenzm**gam)/dz

c       lnuc=ipart(kchg(jamcomp(k(2,i)),5))
c       force(1,i) = -facv(1,lnuc)*(cof1*rxderv+cof2*dxderv)
c       force(2,i) = -facv(1,lnuc)*(cof1*ryderv+cof2*dyderv)
c       force(3,i) = -facv(1,lnuc)*(cof1*rzderv+cof2*dzderv)

        force(1,i)=force(1,i)-(cof1*rxderv+cof2*dxderv)
        force(2,i)=force(2,i)-(cof1*ryderv+cof2*dyderv)
        force(3,i)=force(3,i)-(cof1*rzderv+cof2*dzderv)

  100 continue

      end

c***********************************************************************

      subroutine jambuuci(msel)

c...Purpose: to  calculate the initial value of coulomb.

      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
c...alp, bet: minimum and  maxmum eigenvalue.
      parameter (alp=0.33d0,bet=5.7d0)
      common/jambc1/rmm(2),nrepeat,mmx,mmz,ncskip

      if(msel.eq.0) then
c....Min. and max. eigenvalues.
        gam=(alp*bet)**0.25d0*sqrt((alp+bet)/2.0d0)
        rmm(1)=gam-sqrt(gam**2-alp*bet)
        rmm(2)=gam+sqrt(gam**2-alp*bet)

c...in: number of iteration.
        nrepeat=10
c...ncskip: number of test particle to be counted
        ncskip=3
c...mmx: maxx-mmx; restricted space for coulomb
        mmx=0
c...mmz: maxz-mmz; restricted space for coulomb
        mmz=0

        return
      endif

c     mbx=maxx-mmx
c     mbz=maxz-mmz
      mxx=maxx-mmx-1
      mxz=maxz-mmz-1


      radta=pard(40)    ! radius of target
      radpr=pard(50)    ! radius of projectile
      masspr=mstd(2)    ! projectile mass
      massta=mstd(5)    ! target mass
      mstapr=mstd(3)    ! proton in target
      msprpr=mstd(6)    ! proton in projectile
      zerta=pard(39)    ! z-position of target
      zerpr=pard(49)    ! z-position of projectile
      xerta=pard(37)    ! x-position of target
      xerpr=pard(47)    ! x-position of projectile
      alpe=parc(23)
      rho0=parc(21)
      pi=paru(1)

c...Initial value of cop(ix,iy,iz) from uniform charge distrib.
      do 100 ix=-mxx,mxx
      do 100 iy=-mxx,mxx
      do 100 iz=-mxz,mxz
         if(massta.ne.0) then
         rtag=sqrt((dble(ix)-xerta)**2+dble(iy)**2+
     &             (dble(iz)-zerta)**2)

         if(rtag.lt.radta) then
           cop(ix,iy,iz)=alpe*4d0*pi*rho0
     $     *(radta**2/2d0-rtag**2/6d0)*dble(mstapr)/dble(massta)
         else
           cop(ix,iy,iz)=alpe*dble(mstapr)/rtag
         end if
         end if

         if(masspr.ne.0) then
           rpro=sqrt((dble(ix)-xerpr)**2+dble(iy)**2
     $                                   +(dble(iz)-zerpr)**2)
           if(rpro.lt.radpr) then
             cop(ix,iy,iz)=cop(ix,iy,iz)+
     &          alpe*4d0*pi*rho0*(radpr**2/2d0-rpro**2/6d0)
     &                 *dble(msprpr)/dble(masspr)
           else
             cop(ix,iy,iz)=cop(ix,iy,iz)+alpe*dble(msprpr)/rpro
           end if
         end if
  100 continue

      call jambuucb
      call jambuucl
      nrepeat=mstc(102)

      end

c***********************************************************************

      subroutine jambuucb

c...Purpose: to calculate the coulomb potential on the boundary.

      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      dimension pcop(-maxz:maxz)
      common/jambc1/rmm(2),nrepeat,mmx,mmz,ncskip

      mbx=maxx-mmx
      mbz=maxz-mmz
c     mxx=maxx-mmx-1
c     mxz=maxz-mmz-1
      alpe=parc(23)

      do 101 iix=-mbx,mbx,2*mbx
      do 102 iiy=-mbx,mbx,2
      do 103 iiz=-mbz,mbz,2
        sek=0.0d0

c....Loop over all ensumble.
        do 202 i=1,nv
          if(k(1,i).gt.10) goto 202
          if(pard(1).lt.r(4,i)) goto 202
          kch=jamchge(k(2,i))
          if(mod(k(8,i)-1,ncskip).ne.0.or.kch.eq.0) goto 202
          dis=sqrt((r(1,i)-dble(iix))**2
     &            +(r(2,i)-dble(iiy))**2
     &            +(r(3,i)-dble(iiz))**2)
          if(dis.gt.1.d0) sek=sek+kch/dis/3
  202 continue

      cop(iix,iiy,iiz)=alpe*sek/dble(mstc(5)/ncskip)
      if(iiz.ne.-mbz) then
        cop(iix,iiy,iiz-1)=(cop(iix,iiy,iiz)+prcop)/2.0d0
      end if
      if(iix.ne.-mbx) then
        cop(iix-1,iiy,iiz)=(cop(iix,iiy,iiz)+pcop(iiz))/2.0d0
        if(iiz.ne.-mbz) then
         cop(iix-1,iiy,iiz-1)=(cop(iix,iiy,iiz-1)+pcop(iiz-1))/2.0d0
        end if
      end if
      if(iiz.ne.-mbz) then
        pcop(iiz-1)=cop(iix,iiy,iiz-1)
      end if
      prcop=cop(iix,iiy,iiz)
      pcop(iiz)=cop(iix,iiy,iiz)

  103 continue
  102 continue
  101 continue


      do 111 iiy=-mbx,mbx,2*mbx
      do 112 iix=-mbx,mbx,2
      do 113 iiz=-mbz,mbz,2

        sek=0.0d0
        do 212 i=1,nv
          if(k(1,i).gt.10) goto 212
          if(pard(1).lt.r(4,i)) goto 212
          kch=jamchge(k(2,i))
          if(mod(k(8,i)-1,ncskip).ne.0.or.kch.eq.0) goto 212
          dis=sqrt((r(1,i)-dble(iix))**2
     &            +(r(2,i)-dble(iiy))**2
     &            +(r(3,i)-dble(iiz))**2)
          if(dis.gt.1.d0) sek=sek+kch/dis/3
  212 continue
      cop(iix,iiy,iiz)=alpe*sek/dble(mstc(5)/ncskip)
      if(iiz.ne.-mbz) then
        cop(iix,iiy,iiz-1)=(cop(iix,iiy,iiz)+prcop)/2.0d0
      end if
      if(iix.ne.-mbx) then
        cop(iix-1,iiy,iiz)=(cop(iix,iiy,iiz)+pcop(iiz))/2.0d0
        if(iiz.ne.-mbz) then
         cop(iix-1,iiy,iiz-1)=(cop(iix,iiy,iiz-1)+pcop(iiz-1))/2.0d0
        end if
      end if
      if(iiz.ne.-mbz) then
        pcop(iiz-1)=cop(iix,iiy,iiz-1)
      end if
      prcop=cop(iix,iiy,iiz)
      pcop(iiz)=cop(iix,iiy,iiz)

  113 continue
  112 continue
  111 continue


      do 121 iiz=-mbz,mbz,2*mbz
      do 122 iix=-mbx,mbx,2
      do 123 iiy=-mbx,mbx,2
        sek=0.0d0
        do 222 i=1,nv
          if(k(1,i).gt.10) goto 222
          if(pard(1).lt.r(4,i)) goto 222
          kch=jamchge(k(2,i))
          if(mod(k(8,i)-1,ncskip).ne.0.or.kch.eq.0) goto 222
          dis=sqrt((r(1,i)-dble(iix))**2
     &            +(r(2,i)-dble(iiy))**2
     &            +(r(3,i)-dble(iiz))**2)
          if(dis.gt.1.d0) sek=sek+kch/dis/3
  222 continue

      cop(iix,iiy,iiz)=alpe*sek/dble(mstc(5)/ncskip)
      if(iiy.ne.-mbx) then
        cop(iix,iiy-1,iiz)=(cop(iix,iiy,iiz)+prcop)/2.0d0
      end if
      if(iix.ne.-mbx) then
        cop(iix-1,iiy,iiz)=(cop(iix,iiy,iiz)+pcop(iiy))/2.0d0
        if(iiy.ne.-mbx) then
         cop(iix-1,iiy-1,iiz)=(cop(iix,iiy-1,iiz)+pcop(iiy-1))/2.0d0
        end if
      end if
      if(iiy.ne.-mbx) then
        pcop(iiy-1)=cop(iix,iiy-1,iiz)
      end if
      prcop=cop(iix,iiy,iiz)
      pcop(iiy)=cop(iix,iiy,iiz)

  123 continue
  122 continue
  121 continue

      end

c***********************************************************************

      subroutine jambuucl

c...Purpose: to calculate the coulomb potential from (m) to (m+1)
c...of the iteration.

      include 'jam1.inc'
      include 'jam2.inc'
      include 'jam3.inc'
      dimension w(0:mxv),g(0:mxv),ee(0:mxv)
      common/jambc1/rmm(2),nrepeat,mmx,mmz,ncskip

      mxx=maxx-mmx-1
      mxz=maxz-mmz-1
      mrepeat=2*nrepeat
      alpe=parc(23)

      do 1000 j=1,mrepeat

      ii=1
      jj=j/2*2
      if(jj.eq.j) ii=2

      w(0)=0.0d0
      g(0)=0.0d0
      i=0
      do 100 iy=-mxx,mxx
      do 100 iz=-mxz,mxz
      do 100 ix=-mxx,mxx

         i=i+1
         ee(i)=alpe*4.0d0*paru(1)*rhoz(ix,iy,iz)
     $        +(rmm(ii)-4.0d0)*cop(ix,iy,iz)
     &        + cop(ix  ,iy+1,iz  )
     &        + cop(ix  ,iy-1,iz  )
     &        + cop(ix  ,iy  ,iz+1)
     &        + cop(ix  ,iy  ,iz-1)
         if(ix.eq.-mxx) then
             ee(i)=ee(i)+cop(ix-1,iy,iz)
             ai=0.0d0
         else
             ai=-1.0d0
         end if
         if(ix.eq. mxx) then
             ee(i)=ee(i)+cop(ix+1,iy,iz)
             ci=0.0d0
         else
             ci=-1.0d0
         end if
         g(i)=(ee(i)-ai*g(i-1))/(2.0d0+rmm(ii)-ai*w(i-1))
         w(i)=ci/(2.0d0+rmm(ii)-ai*w(i-1))
  100 continue

         i=i+1
         coppre=0.0d0
      do 200 iy=mxx,-mxx,-1
      do 200 iz=mxz,-mxz,-1
      do 200 ix=mxx,-mxx,-1
         i=i-1
         cop(ix,iy,iz)=g(i)-w(i)*coppre
         coppre=cop(ix,iy,iz)
  200 continue

c-----------------------------------------------------------------------
         i=i-1
      do 101 iz=-mxz,mxz
      do 101 ix=-mxx,mxx
      do 101 iy=-mxx,mxx
         i=i+1
         ee(i)=alpe*4.0d0*paru(1)*rhoz(ix,iy,iz)
     $        +(rmm(ii)-4.0d0)*cop(ix,iy,iz)
     &        + cop(ix+1,iy  ,iz  )
     &        + cop(ix-1,iy  ,iz  )
     &        + cop(ix  ,iy  ,iz+1)
     &        + cop(ix  ,iy  ,iz-1)
         if(iy.eq.-mxx) then
             ee(i)=ee(i)+cop(ix ,iy-1,iz )
             ai=0.0d0
         else
             ai=-1.0d0
         end if
         if(iy.eq. mxx) then
             ee(i)=ee(i)+cop(ix ,iy-1,iz )
             ci=0.0d0
         else
             ci=-1.0d0
         end if
         g(i)=(ee(i)-ai*g(i-1))/(2.0d0+rmm(ii)-ai*w(i-1))
         w(i)=ci/(2.0d0+rmm(ii)-ai*w(i-1))
  101 continue
         i=i+1
         coppre=0.0d0
      do 201 iz=mxz,-mxz,-1
      do 201 ix=mxx,-mxx,-1
      do 201 iy=mxx,-mxx,-1
         i=i-1
         cop(ix,iy,iz)=g(i)-w(i)*coppre
         coppre=cop(ix,iy,iz)
  201 continue


         i=i-1
      do 102 ix=-mxx,mxx
      do 102 iy=-mxx,mxx
      do 102 iz=-mxz,mxz
         i=i+1
         ee(i)=alpe*4.0d0*paru(1)*rhoz(ix,iy,iz)
     $        +(rmm(ii)-4.0d0)*cop(ix,iy,iz)
     &        + cop(ix+1,iy  ,iz)
     &        + cop(ix-1,iy  ,iz)
     &        + cop(ix  ,iy+1,iz)
     &        + cop(ix  ,iy-1,iz)
         if(iz.eq.-mxz) then
             ee(i)=ee(i)+cop(ix,iy,iz-1)
             ai=0.0d0
         else
             ai=-1.0d0
         end if
         if(iz.eq. mxz) then
             ee(i)=ee(i)+cop(ix,iy,iz+1)
             ci=0.0d0
         else
             ci=-1.0d0
         end if
         g(i)=(ee(i)-ai*g(i-1))/(2.0d0+rmm(ii)-ai*w(i-1))
         w(i)=ci/(2.0d0+rmm(ii)-ai*w(i-1))
  102 continue
         i=i+1
         coppre=0.0d0
      do 202 ix=mxx,-mxx,-1
      do 202 iy=mxx,-mxx,-1
      do 202 iz=mxz,-mxz,-1
         i=i-1
         cop(ix,iy,iz)=g(i)-w(i)*coppre
         coppre=cop(ix,iy,iz)
  202 continue

 1000 continue  ! iteration loop closed.

      end
