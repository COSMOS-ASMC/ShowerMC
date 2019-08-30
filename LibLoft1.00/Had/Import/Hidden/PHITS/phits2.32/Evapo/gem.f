************************************************************************
*                                                                      *
      subroutine gemset
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              initialize GEM for PHITS                                *
*                                                                      *
*        parameters :                                                  *
*                                                                      *
*        alev    = 0.0  ; The GCCI level density parameter is used.    *
*                = 1.0  ; The level density parameter is given by a=A/8*
*                > 1.0  ; The level density parameter.                 *
*                         is given by a = A / alev.                    *
*                                                                      *
*        rcal    = 0.0  ; The Dostrovsky and Matsuse parameter set     *
*                         are used for the inverse reaction            *
*                         cross sections.                              *
*                = 10.0 ; The simple parameter set is used for         *
*                         the inverse reaction cross section           *
*                         r0 is set to 1.5.                            *
*                < 10.0 ; The simple parameter set is used for         *
*                         the inverse reaction cross section           *
*                         r0 is set to rcal.                           *
*                                                                      *
*        ifis    ne 0   ; Original parameter in the Atchison model     *
*                         is used.                                     *
*                = 0    ; New parameter set is used.                   *
*                                                                      *
*        nimax   : maximum number of ejectiles                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      common /options/alev, rcal, ifis
      common /exiejn/ nimax,ndmax,nlvl(70)

*-----------------------------------------------------------------------
*     read parameters
*-----------------------------------------------------------------------

            alev = 0.0d0
            rcal = 0.0d0
            ifis = 0

*-----------------------------------------------------------------------

            call setup

            nimax = ndmax
c           nimax = 6

*-----------------------------------------------------------------------

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine gemexec(iz,in,ex,px,py,pz,pt,et,rm,wt,ierr)
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to execute statistical particle decay / fission.        *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*        iz, in     : proton and neutron number of mather              *
*        ex         : excitation energy of mather (MeV)                *
*        px,py,pz   : momentum vector of mather (GeV)                  *
*        pt         : absolute value of momentum of mather (GeV)       *
*        et         : energy of mather (sqrt(p**2+m**2) GeV)           *
*        rm         : rest mass of mather (GeV)                        *
*        wt         : weight change                                    *
*        ierr       : error flag                                       *
*                                                                      *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'

*-----------------------------------------------------------------------

      parameter ( amu = 931.494d0 )

*-----------------------------------------------------------------------

      common /clusts/ nclust, kclust(3,nnn)
      common /clustu/ lclust(0:7,nnn), sclust(0:12,nnn)
      common /clustv/ kdecay(4)
      common /emode/ ge1, ge2, iemode, nevhin, nlowlv, nefiss, ntwidt

*-----------------------------------------------------------------------

      common /mom/ erec,bet0(3),ekin,bet1(3)
      dimension ir(2), ir1(2)

      common /sim/nocas
      common /fiss3/ zfis(2),afis(2),ufis(2),er(2),betf(2,3)

      logical fisinh

      common /gemint/ initgm
      common /gamint/ initga

*-----------------------------------------------------------------------
*     initialization
*-----------------------------------------------------------------------

         if( initgm .eq. 0 ) then

            initgm = initgm + 1

            call gemset

         end if

*-----------------------------------------------------------------------
*     initialization for gamlib
*-----------------------------------------------------------------------

         if( initga .eq. 0 ) then

            initga = initga + 1

c/////            call binin2(ierr)   ! since igamma=0, we don't consider gamma
                                      ! from deexcited A

c////               if( ierr .ne. 0 ) call parastop( 789 )

         end if

*-----------------------------------------------------------------------
*        initial nucleus and values
*-----------------------------------------------------------------------

               ierr = 0

               erec = ( et - rm ) * 1000.0

            if( pt .gt. 0.0d0 ) then

               bet0(1) = px / pt
               bet0(2) = py / pt
               bet0(3) = pz / pt

            else

               bet0(1) = 0.0
               bet0(2) = 0.0
               bet0(3) = 1.0

            end if

               iz0 = iz
               in0 = in
               ia0 = in + iz

               e1  = ex
               a0  = dble(ia0)
               z0  = dble(iz0)

               iflag = 0
               nocas = 0

*-----------------------------------------------------------------------
*        initialize fission fragment
*-----------------------------------------------------------------------

          fisinh = .false.

          do k = 1, 2

            afis(k) = 0.0d0
            zfis(k) = 0.0d0
            ufis(k) = 0.0d0
            er(k)   = 0.0d0

            do kk = 1, 3
              betf(k,kk) = 0.d0
            end do

          end do

*-----------------------------------------------------------------------
*     Start decay.
*        zero set for booking
*-----------------------------------------------------------------------

               nclust = 0
               i1     = 0

               ifssev = 0

            do i = 1, 4

               kdecay(i) = 0

            end do

*-----------------------------------------------------------------------

            imstp = 0

 2000    continue

            imstp = imstp + 1

            if( i1 + 2 .gt. nnn ) goto 5000

*-----------------------------------------------------------------------
*        for event generator mode 
*-----------------------------------------------------------------------

               ntwidt = nevhin

            if( nevhin .ge. 2 .and. imstp .ge. nevhin ) then

               ntwidt = 1

            end if

*-----------------------------------------------------------------------
*        start evaporation calculation
*
*        iflag=1 : no more emission
*        iflag=2 : fission
*        iflag=3 : emission occured
*-----------------------------------------------------------------------

            call stdcay(a0,z0,e1,iaf,izf,iflag,fisinh)

            if( .not. fisinh ) goto 40

               ifssev = ifssev + 1

*-----------------------------------------------------------------------
*     fission fragment and ejectile from fission fragment
*-----------------------------------------------------------------------

         do 510 k = 1, 2

               a0   = afis(k)
               z0   = zfis(k)
               e1   = ufis(k)
               erec = er(k)

               do kk = 1, 3
                  bet0(kk) = betf(k,kk)
               end do

 3000    continue

         if( i1 + 2 .gt. nnn ) goto 5000

               e1old = e1
               aold  = a0
               zold  = z0

               call stdcay(a0,z0,e1,iaf,izf,iflag,fisinh)

            if( iflag .eq. 3 ) then

               i1 = i1 + 1
               erex = 0.0
               call gemout(i1,izf,iaf,ekin,bet1,erex,wt)


            else if( iflag .eq. 1 )then

               izr = nint(z0)
               iar = nint(a0)

               i1 = i1 + 1
               call gemout(i1,izr,iar,erec,bet0,e1,wt)

               goto 510

            end if

               goto 3000

 510    continue

               goto 6000

*-----------------------------------------------------------------------
*     ejectile and residual nucleus without fission
*-----------------------------------------------------------------------

 40      if( iflag .eq. 3 ) then

               i1 = i1 + 1
               erex = 0.0
               call gemout(i1,izf,iaf,ekin,bet1,erex,wt)

         else if( iflag .eq. 1 )then

               izr = nint(z0)
               iar = nint(a0)

               i1 = i1 + 1
               call gemout(i1,izr,iar,erec,bet0,e1,wt)

               goto 6000

         end if

*-----------------------------------------------------------------------

      goto 2000

*-----------------------------------------------------------------------
*     end of evaporation
*-----------------------------------------------------------------------

 6000    continue

            if( ifssev .gt. 0 )  kdecay(4) = 1

         return

*-----------------------------------------------------------------------
*           check dimension
*-----------------------------------------------------------------------

 5000    continue

               write(*,*) ' **** Error at gemexec, too many products'
               write(*,*) ' ========================================='
               write(*,*) ' nnn  = ', nnn

               ierr = 1
               call parastop( 999 )

*-----------------------------------------------------------------------

      return
      end  subroutine


************************************************************************
*                                                                      *
      subroutine gemout(i1,iz,ia,ekin0,bet,ex,wt)
*                                                                      *
*                                                                      *
*       booking of out going particles from GEM
*                                                                      *
*     input:                                                           *
*                                                                      *
*        i1         : number of ejectile                               *
*        iz, ia     : proton and mass number of ejectile               *
*        ekin0      : kinetic energy of ejectile (MeV)                 *
*        bet        : momentum unit vector                             *
*        ex         : excitation energy of ejectile (MeV)              *
*        wt         : weight change                                    *
*                                                                      *
*     output:                                                          *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclust   : total number of out going particles and nuclei     *
*                                                                      *
*        kclust(3,nclust)                                              *
*                                                                      *
*                   kclust(1,i) = 101 : final output                   *
*                   kclust(2,i) = 0                                    *
*                   kclust(3,i) = 0                                    *
*                                                                      *
*        lclust(i,nclust)                                              *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3,                                                *
*                  = 4,                                                *
*                  = 5, charge                                         *
*                  = 6,                                                *
*                  = 7,                                                *
*                                                                      *
*        sclust(i,nclust)                                              *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, px (GeV/c)                                     *
*                  = 2, py (GeV/c)                                     *
*                  = 3, pz (GeV/c)                                     *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight change                                  *
*                  = 9, delay time                                     *
*                  = 10, x-displace                                    *
*                  = 11, y-displace                                    *
*                  = 12, z-displace                                    *
*                                                                      *
*        kdecay(4) = 0 : no fission                                    *
*                  = 1 : with fission                                  *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( amu = 931.494d0 )

*-----------------------------------------------------------------------
*     common for output of evaporation and fission
*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clusts/ nclust, kclust(3,nnn)
      common /clustu/ lclust(0:7,nnn), sclust(0:12,nnn)
      common /clustv/ kdecay(4)

      dimension bet(3)

*-----------------------------------------------------------------------
*     number of ejectiles
*-----------------------------------------------------------------------

            nclust = i1

*-----------------------------------------------------------------------

            ekin = ekin0 / 1000.0
            rmsp = ( amu * dble( ia ) + energm(iz,ia) ) / 1000.0
            pabs = sqrt( ekin**2 + 2.d0 * rmsp * ekin )

*-----------------------------------------------------------------------

               kclust(1,nclust)  = 101
               kclust(2,nclust)  = 0
               kclust(3,nclust)  = 0

               lclust(0,nclust)  = 0
               lclust(1,nclust)  = iz
               lclust(2,nclust)  = ia - iz
               lclust(3,nclust)  = 0
               lclust(4,nclust)  = 0
               lclust(5,nclust)  = iz
               lclust(6,nclust)  = 0
               lclust(7,nclust)  = 0

               sclust(0,nclust)  = 0.0
               sclust(1,nclust)  = pabs * bet(1)
               sclust(2,nclust)  = pabs * bet(2)
               sclust(3,nclust)  = pabs * bet(3)
               sclust(4,nclust)  = ekin + rmsp
               sclust(5,nclust)  = rmsp
               sclust(6,nclust)  = ex
               sclust(7,nclust)  = ekin * 1000.
               sclust(8,nclust)  = wt
               sclust(9,nclust)  = 0.0
               sclust(10,nclust) = 0.0d0
               sclust(11,nclust) = 0.0d0
               sclust(12,nclust) = 0.0d0

*-----------------------------------------------------------------------

      return
      end subroutine



************************************************************************
*                                                                      *
      subroutine setup
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              SETUP Nuclear Data                                      *
*                 Excess Mass Table                                    *
*                 Nuclear levels of the ejectiles                      *
*                 Shell effect data for Myer-Swaitecki fission barrier *
*                 calculation                                          *
*                                                                      *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)

      parameter ( nmastb = 9796 )
      parameter ( nshltb = 9107 )

*-----------------------------------------------------------------------

      common /exims/ wapsm(0:150,0:250)
      common /exiejn/ nimax,ndmax,nlvl(70)
      common /exieje/ exm(70,200),spin(70,200),width(70,200)
      common /ejectl/ omega(70),ifa(70),ifz(70)
      common /fiss4/ shell(150,250)

      common /gemms/ niz(10000), nia(10000), bie(10000)
      common /gemsh/ nsz(10000), nsa(10000), shin(10000)

      data wapsm / 37901*0.0d0 /
      data shell / 37500*0.0d0 /

      character cname*100
      character nname*4

*-----------------------------------------------------------------------
*     Read excess mass table
*-----------------------------------------------------------------------

            do j = 1, nmastb

               iz = niz(j)
               ia = nia(j)
               wapsm(iz,ia-iz) = bie(j)

            end do

*-----------------------------------------------------------------------
*     Read data for the excitation state of the ejectiles
*-----------------------------------------------------------------------

            do i = 1, 6

                  nlvl(i) = 0

            do j = 1, 200
                   exm(i,j) = 0.d0
                  spin(i,j) = 0.d0
                 width(i,j) = 0.d0
            end do
            end do

            do i = 7, ndmax
            do j = 1, 200

               if( j .gt. nlvl(i) ) then

                   exm(i,j) = 0.d0
                  spin(i,j) = 0.d0
                 width(i,j) = 0.d0

               end if

            end do
            end do

*-----------------------------------------------------------------------
*     Read shell effect data for fission barrier calculations
*-----------------------------------------------------------------------

            do j = 1, nshltb

               iz = nsz(j)
               ia = nsa(j)
               shell(iz,ia-iz) = shin(j)

            end do

*-----------------------------------------------------------------------

      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine stdcay(a,z,u,iaf,izf,iflag,fisinh)
C/////////////////////////////////////////////////////////////////////
C <STDCAY>
C  Calculate evaporation process and fission process
C    - determine if emission occurs
C    - determine if fission occurs
C        -if fission , call FIS
C        -if no fission, determine kinetic energy, etc. for emittor
C                        and caluculate recoil energy, etc.
C=====================================================================
C <Subroutine>
C   gamma    :  Decay width calculation
C    fis     :  Fission calculation
C selectE    :  Select kinetic energy in the CM system
C--------------------------------------------------------------------
C <Function>
C  fprob     : Calculate fission probabirity
C=====================================================================
C <Variables>
C     a    : mass of parent nuclei -> residual nuclei   (IN and OUT)
C     z    : charge of parent nuclei    -> residual     (IN and OUT)
C     u    : excited energy of parent nuclei-> residual (IN and OUT)
C   iflag  : 1 = no more emission                       (OUT)
C          : 2 = fission
C          : 3 = emission occur
C  fisinh  : true=fission occur, false=no fission       (OUT)
c   erec   : recoil energy in the lab system            (IN and OUT)
c   ekin   : kinetic energy in the lab system           (OUT)
c   bet0   : unit vector of recoil momentum             (IN and OUT)
c   bet1   : unit vector of momentum of emittor         (IN and OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)

      parameter (amu=931.494d0)
      parameter (pi=3.1415926535898d0)

      common /mom/ erec,bet0(3),ekin,bet1(3)

      common /exiejn/ nimax,ndmax,nlvl(70)

      common /std1/ r(70),s(70),sigma,rr(70)
      common /ejectl/ omega(70),ifa(70),ifz(70)
      common /emitr/ gj(70),q(70),V(70),delta(70),smalla(70)

      common /fiss3/ zfis(2),afis(2),ufis(2),er(2),betf(2,3)

*-----------------------------------------------------------------------
*     common for gamlib
*-----------------------------------------------------------------------

      common /qparm/ ielas,icasc,iqstep,lvlopt,igamma
      common /emode/ ge1, ge2, iemode, nevhin, nlowlv, nefiss, ntwidt

      parameter (len1 =41600)
      parameter (len2 = 1600)

      common/ big1   /wtnone,izao(len2),ipt(len2),iref(105)
      common/ big2   /elo(len1),isp(len1)

*-----------------------------------------------------------------------

      logical fisinh

C  Initialization
      iaf=0
      izf=0
      ia=nint(a)
      iz=nint(z)

C   Evaporation calculation starts
      uran=rn(0)

C  Calculate decay width

      call gammag(a,z,u,uran)

      if (sigma.le.0.0d0) then
       iflag=1  ! no more emission
       return
      endif

C... Select ejectile
      uran=uran*sigma
      sum=0.d0

c     do 20 j=1,70
      do 20 j = 1, nimax

        k=j
        sum=sum+r(j)
        if (sum.ge.uran) go to 30

 20   continue

c     write(*,*)'ERROR: Ejectile is not selected in subroutime STDCAY:
c    &     sum r(j) < uran ', sum, uran
c     call parastop( 999 )
cKN
       iflag=1
       u = 0.0d0
       return
cKN

 30   continue
      jemiss=k

*-----------------------------------------------------------------------
C  Fission calculation :
C  fission only occurs one time
C  (fission never occur to post-fission nucleus)

      if (jemiss.eq.1.and.z.gt.70.d0.and..not.fisinh) then

       probf=fprob(z,a,u)

*-----------------------------------------------------------------------
cKN 2010/01/26

       if( ( probf .gt. rn(0) .and. nefiss .eq. 0 ) .or. 
     &     ( probf .gt. 0.d0  .and. nefiss .gt. 0 ) ) then

cKN
*-----------------------------------------------------------------------

        fisinh=.true.
        call fis(a,z,u,fisinh)
        fisinh=.true.
        iflag=2                 !: fission
        return
       endif

      endif

*-----------------------------------------------------------------------
C...set A & Z of the ejectile

      iaf=ifa(jemiss)
      izf=ifz(jemiss)
      iflag=3   ! Emission occured

*-----------------------------------------------------------------------
C...select kinetic energy in the CM system
*-----------------------------------------------------------------------

      pran = rn(0)

      pekin = pran * r(jemiss) * rr(jemiss)

        call selectE(jemiss,a,z,u,pekin,ekin)

        if(ekin.le.0.d0) then
c        write(*,*)'ERROR: Kinetic energy can not be determined in
c    &        subroutine selectE'
c        call parastop( 999 )
cKN
       iflag=1
       u = 0.0d0
       return
cKN
        endif

*-----------------------------------------------------------------------
C... Velocity of Parent nucleus in the LAB system: vres*bet0

      am=a*amu+energm(iz,ia)
      pres=sqrt(erec**2+2.d0*am*erec)
      vres=pres/(Erec+am)

*-----------------------------------------------------------------------

C... residual nucleus, A, Z
      a=a-dble(iaf)
      z=z-dble(izf)

*-----------------------------------------------------------------------
*     check level near ground state
*     and adjust excitation energy and exit energy
*-----------------------------------------------------------------------

      if( nlowlv .eq. 1 ) then

                  e = u - q(jemiss) - ekin

         if( e .gt. 0.0d0 ) then

                  ja  = nint( a )
                  jz  = nint( z )

                  jza = 1000*jz+ja

*-----------------------------------------------------------------------
*           special for 10B + n -> 11B -> alpha + 7Li
*           first excited state / ground = 14.9
*-----------------------------------------------------------------------

            if( ja+iaf .eq. 11 .and. jz+izf .eq. 5 .and.
     &          jza .ne. 3007 .and. u - q(jemiss) .lt. 4.63 ) then

                    iflag=1
                    u = 0.0d0
                    return

            else
     &      if( ja+iaf .eq. 11 .and. jz+izf .eq. 5 .and.
     &          jza .eq. 3007 .and. u - q(jemiss) .lt. 4.63 ) then

               if( 15.9 * rn(0) .gt. 1.d0 ) then

                  e = 0.47761d0

               else

                  e = 0.0d0

               end if

*-----------------------------------------------------------------------
*           special for 3He + n -> 4He -> 3H + 1H
*           all to ground
*-----------------------------------------------------------------------

            else if( ja+iaf .eq. 4 .and. jz+izf .eq. 2 .and.
     &               jza .eq. 1003 .and. u - q(jemiss) .lt. 1.00 ) then

                 e = 0.0d0

*-----------------------------------------------------------------------

            else

                  k1 = iref(jz)
                  k2 = iref(jz+1)-1

               do k = k1, k2

                  if( jza .eq. izao(k) ) then

                     ik  = k
                     nup = ipt(ik+1)-1
                     nlo = ipt(ik)

                     do i = nlo + 1, nup

                        if( e .le. elo(i) ) then

                              e = elo(i-1)

                           goto 10

                        end if

                     end do

                  end if

               end do

   10          continue

            end if

         end if

                  ekin = u - q(jemiss) - e

                  if(ekin.le.0.d0) then
                     iflag=1
                     u = 0.0d0
                     return
                  end if

      end if

*-----------------------------------------------------------------------
C... excitation energy

      u = u - q(jemiss) - ekin

*-----------------------------------------------------------------------

C... Momentum in the CM system : pcmx, pcmy, pcmz
      amf=iaf*amu+energm(izf,iaf)
      amr=(ia-iaf)*amu+energm(iz-izf,ia-iaf)
      redm=amf*amr/am

                  if(redm.le.0.d0) then
                     iflag=1
                     u = 0.0d0
                     return
                  end if

      pcm=sqrt(ekin**2+2.d0*redm*ekin)
      th  = acos( 2.0d0 * rn(0) - 1.0d0 )
      ph  = 2.0d0 * pi  * rn(0)

      pcmx = pcm* dsin(th) * dcos(ph)
      pcmy = pcm* dsin(th) * dsin(ph)
      pcmz = pcm* dcos(th)

*-----------------------------------------------------------------------

C... Boost emjectile momentum to Lab system : p1x,p1y,p1z
      gam=sqrt(1.d0-vres**2)
      gam=1.d0/gam
      pv=pcmx*bet0(1)+pcmy*bet0(2)+pcmz*bet0(3)
      pv=pv*vres

      e1cm=sqrt(pcm**2+amf**2)
      tr=gam*(e1cm+gam*pv/(gam+1))

      p1x=pcmx+vres*bet0(1)*tr
      p1y=pcmy+vres*bet0(2)*tr
      p1z=pcmz+vres*bet0(3)*tr

C... Boost residual momentum to Lab system: p2x,p2y,p2z
      e2cm=sqrt(pcm**2+amr**2)
      tr=gam*(e2cm+gam*(pv)/(gam+1))

      p2x=-pcmx+vres*bet0(1)*tr
      p2y=-pcmy+vres*bet0(2)*tr
      p2z=-pcmz+vres*bet0(3)*tr

C... Kinetic Energy and velocity of Residual in Lab system: erec, bet0
      pr2=p2x**2+p2y**2+p2z**2
      erec= sqrt(amr**2+ pr2)-amr

      bet0(1)=p2x/sqrt(pr2)
      bet0(2)=p2y/sqrt(pr2)
      bet0(3)=p2z/sqrt(pr2)

C... Kinetic Energy and velocity of Emittor in Lab system: ekin, bet1
      pe2=p1x**2+p1y**2+p1z**2
      ekin= sqrt(amf**2+ pe2)-amf

      bet1(1)=p1x/sqrt(pe2)
      bet1(2)=p1y/sqrt(pe2)
      bet1(3)=p1z/sqrt(pe2)

      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine gammag(a,z,u,uran)
C/////////////////////////////////////////////////////////////////////
C GAMMA
C  calculate decay width for each particle emission
C=====================================================================
C <Subroutine>
C   eye    : Calculate I0,I1,I2,I3
C--------------------------------------------------------------------
C <Function>
C  dost    : Calculate kp, cp, or k_alpha
C  paire   : Calculate pairing energy
C=====================================================================
C <variables>
C     a    :   residual mass before emission                     (IN)
C     z    :   charge number of res nuclei before emission       (IN)
C     u    :   excitation energy of nuclei before emission       (IN)
C   uran   :   random number                                     (IN)
C     r    :   decay width                                      (OUT)
C     rr   :   decay width enchancement factor by excited state
C                            particle emission                  (OUT)
C   ifa    :   mass of emittor                                  (OUT)
C   ifz    :   charge of emittor                                (OUT)
C  omega   :   spin of emittor                                  (OUT)
C   gj     :   gj in eq.(##)                                    (OUT)
C   q      :   Q-value                                          (OUT)
C  delta   :   pairing energy                                   (OUT)
C  smalla  :   level density parameter                          (OUT)
C   V      :   Coulomb barrier                                  (OUT)
C   exm    :   mass of an excited state [MeV]                   (OUT)
C  spin    :   spin of an excited state                         (OUT)
C  width   :   lifetime of an excited state [MeV]               (OUT)
C   gamn   :   Decay width for neutron emission                 (OUT)
C   an     :   Level density parameter for neutron emission     (OUT)
C beta,alp :  Inverse cross section parameters for neutron emission(OUT)
C/////////////////////////////////////////////////////////////////////
      implicit doubleprecision(a-h,o-z)

      common /ejectl/omega(70),ifa(70),ifz(70)
      common /emitr/ gj(70),q(70),V(70),delta(70),smalla(70)
      common /exiejn/ nimax,ndmax,nlvl(70)
      common /exieje/ exm(70,200),spin(70,200),width(70,200)
      common /std1/ r(70),s(70),sigma,rr(70)
      common /fiss/ beta,alp
      common /options/alev, rcal, ifis

      common /sim/nocas
cKN
      common /emode/ ge1, ge2, iemode, nevhin, nlowlv, nefiss, ntwidt
cKN
      real*8 hbarc
      parameter (hbarc=197.327053d0)
      parameter (amu=931.494d0)
      parameter (pi=3.1415926535898d0)

      dimension couk(70),couc(70)
      data couk /70*1.d0/
      data couc /70*0.d0/
      save couk,couc !FURUTA

*-----------------------------------------------------------------------
*        Parent nucleus
*-----------------------------------------------------------------------

               ia = nint(a)
               iz = nint(z)
               q1 = energm(iz,ia)

*-----------------------------------------------------------------------

      if(q1.eq.1.d10) then

c      write(*,*)'ERROR in GEM: Mass calculation error occurred',
c    & 'in function ENEGY for the nucleus with iz= ',iz, ' ia= ',ia
c        call parastop( 999 )
cKN
         sigma = -1.0
         return
cKN
      endif

*-----------------------------------------------------------------------
*        Set cj and kj for Dostrovsky's parameter set
*-----------------------------------------------------------------------

            if( rcal .eq. 0.0 ) then

               couk(2) = dostg(1,z-ifz(2))

               couk(3) = couk(2) + 0.06
               couk(4) = couk(2) + 0.12

               couk(6) = dostg(2,z-ifz(6))
               couk(5) = couk(6) - 0.06


               couc(2) = dostg(3,z-ifz(2))

               couc(3) = couc(2) / 2.d0
               couc(4) = couc(2) / 3.d0

            end if

*-----------------------------------------------------------------------
*     Start calcualting Gamma for each particle emission
*-----------------------------------------------------------------------

      do 1 j = 1, nimax

*-----------------------------------------------------------------------
*        initialization
*-----------------------------------------------------------------------

             s(j)  = 0.0d0
             r(j)  = 0.0d0
             rr(j) = 1.0d0
             gj(j) = 0.0d0

             if( ifa(j) .eq. 0 ) goto 1

*-----------------------------------------------------------------------
*        special calculation for K.NIITA
*        without neutron width
*-----------------------------------------------------------------------
cKN
         if( ntwidt .eq. 1 .and. j .eq. 1 ) goto 1
         if( ntwidt .ge. 2 .and. j .gt. 1 ) goto 1

*-----------------------------------------------------------------------
*        Reduce calculation time by KN
*-----------------------------------------------------------------------

            if( j .gt. 6 ) then

               if( ia .gt. 40 .and.
     &             uran .lt. 0.95 ) goto 1

               if( ia .gt. 30 .and. ia .le. 40 .and.
     &             uran .lt. 0.93 ) goto 1

               if( ia .gt. 20 .and. ia .le. 30 .and.
     &             uran .lt. 0.7 ) goto 1

            end if

*-----------------------------------------------------------------------
*        daughter nucleus mass and charge
*-----------------------------------------------------------------------

            iaa = ia - ifa(j)
            aa  = dble(iaa)
            izz = iz - ifz(j)
            zz  = dble(izz)
            nn  = iaa - izz

*-----------------------------------------------------------------------
*        Check of the residual nuclei after the emission
*        and avoid double counting modified by KN
*-----------------------------------------------------------------------
cKN
            if( ntwidt .eq. 1 .and. iaa .eq. 1 .and. nn .eq. 1 ) goto 1

*-----------------------------------------------------------------------

            if( iaa .le. 0 .or.
     &          izz .lt. 0 .or.
     &          iaa .lt. izz ) go to 1

            if( iaa .lt. ifa(j) .or.
     &          izz .lt. ifz(j) ) then

               do k = 1, j - 1

                  if( iaa .eq. ifa(k) .and.
     &                izz .eq. ifz(k) ) goto 1

               end do

             end if

*-----------------------------------------------------------------------
*        Q-value
*-----------------------------------------------------------------------

            q2   = energm(izz,iaa)
            q(j) = q2 - q1 + energm(ifz(j),ifa(j))

*-----------------------------------------------------------------------

          if(q2.eq.1.d10) then

c          write(*,*)'ERROR in GEM: Mass calculation error occurred',
c    &   ' in function ENEGY for the nucleus with iz= ',izz,' ia= ',iaa
c           call parastop( 999 )
cKN
         sigma = -1.0
         return
cKN
          endif

*-----------------------------------------------------------------------
*        Coulomb potential
*-----------------------------------------------------------------------

            V(j) = vcoul(zz,aa,ifz(j),ifa(j),couk(j),j)

*-----------------------------------------------------------------------
*        dependence of excitation energy : temp by KN
*-----------------------------------------------------------------------

            if( j .gt. 1 ) then

               V(j) = V(j) / ( 1.d0 + 0.005d0 * u / dble(ifa(j)) )

            end if

*-----------------------------------------------------------------------
*        Set Coulomb potential to zero if Q-value < 0, 
*        for light residual e.g. 8Be, 9B
*-----------------------------------------------------------------------

            if( q(j) .le. 0 .and. aa .le. 20 ) V(j) = 0.d0

*-----------------------------------------------------------------------
*        Reduce Coulomb potential if neutron width is zero 
*-----------------------------------------------------------------------
cKN

            if( ntwidt .eq. 1 ) then
               if( ia .le. 24 ) v(j) = 0.5 * v(j)
               if( ia .le. 13 ) v(j) = 0.d0
            end if

*-----------------------------------------------------------------------
*        Paring energy
*-----------------------------------------------------------------------

            delta(j) = paire(izz,nn)

*-----------------------------------------------------------------------
*        Check whether the emission is enegetically possible or not
*-----------------------------------------------------------------------

            if( u - q(j) - V(j) .le. 0.d0 ) goto 1

*-----------------------------------------------------------------------
*        alpha and beta for neutrons: the precise parameter set
*-----------------------------------------------------------------------

            if( rcal .eq. 0.0 ) then

               alp  = 0.76d0 + 1.93d0 / aa**.333333333d0
               beta = ( 1.66d0 / aa**.666666666667d0 -5.d-2 ) / alp
               beta = max(beta,0.d0)

            else

               alp  = 1.d0
               beta = 0.d0

            end if

               bett = -V(j)
               if( j .eq. 1 ) bett = beta

*-----------------------------------------------------------------------
*        level density parameter and gamma
*-----------------------------------------------------------------------

            smalla(j) = getag(u-q(j)-V(j),izz,nn,isdum)

            call eye(j,aa,izz,nn,u,q(j),V(j),delta(j),smalla(j),bett,
     &               r(j),s(j))

*-----------------------------------------------------------------------
*        (2Sj+1)mj * alpha
*-----------------------------------------------------------------------

            gj(j) = ( 2.d0 * omega(j) + 1.d0 ) * dble(ifa(j))

         if( j .eq. 1 ) then

            gj(j) = gj(j) * alp

         else

            gj(j) = gj(j) * ( 1.d0 + couc(j) )

         endif

*-----------------------------------------------------------------------
*        geometric cross section... sigma R
*-----------------------------------------------------------------------

            rmass = rb(aa,ifa(j),j)
            gj(j) = rmass * rmass * gj(j)

*-----------------------------------------------------------------------
*        decay width
*-----------------------------------------------------------------------

            r(j) = gj(j) * r(j)
            r(j) = max(r(j),0.d0)

*-----------------------------------------------------------------------
*     decay width calculation from excited states
*-----------------------------------------------------------------------

            if( j .le. 6 .or.
     &          r(j) .eq. 0 .or.
     &          nlvl(j) .eq. 0 ) goto 1

               rrq = r(j)

*-----------------------------------------------------------------------
*        paring eneryg for a parent nucleus
*-----------------------------------------------------------------------

               del    = paire(iz,ia-iz)
               aparnt = getag(u,iz,ia-iz,isdum)

*-----------------------------------------------------------------------
*        level density of mother : rho_i
*-----------------------------------------------------------------------

               rhop = rho(ia,iz,u,del,aparnt)

*-----------------------------------------------------------------------
*        sum up for excited states
*-----------------------------------------------------------------------

         do 10 i = 1, nlvl(j)

*-----------------------------------------------------------------------
*           Q-value for an excited state
*-----------------------------------------------------------------------

               qq = q(j) + exm(j,i)

               if( u - qq - V(j) .le. 0.d0 ) goto 10

*-----------------------------------------------------------------------
*           level denisty for an excited state
*-----------------------------------------------------------------------

               aq = getag(u-qq-V(j),izz,nn,isdum)

*-----------------------------------------------------------------------
*           integral part of gamma for an excited state
*-----------------------------------------------------------------------

               call eye(j,aa,izz,nn,u,qq,V(j),delta(j),aq,bett,rq,sq)

*-----------------------------------------------------------------------
*           (2Sj+1)mj * alpha * sigma_R for an excited state
*-----------------------------------------------------------------------

               gq = ( 2.d0 * spin(j,i) + 1.d0 ) * dble(ifa(j))
               gq = rmass * rmass * gq

*-----------------------------------------------------------------------
*           gamma for an excited state [MeV]
*           reject an excited state
*           if the decay width [MeV] < level widht [MeV]
*-----------------------------------------------------------------------

               rrqg = ( gq * rq / rhop ) * amu / hbarc**2 / pi

               if( width(j,i) .ge. rrqg ) goto 20

*-----------------------------------------------------------------------

               rrq = rrq + gq * rq

  10     continue

*-----------------------------------------------------------------------
*           rr(j) : gamma(1-6) / total gamma
*-----------------------------------------------------------------------

  20     continue

               rr(j) = r(j) / rrq
               r(j)  = rrq

*-----------------------------------------------------------------------

    1 continue

*-----------------------------------------------------------------------
*        total gamma
*-----------------------------------------------------------------------

  900    continue

            sigma = 0.d0

         do j = 1, nimax

            sigma = sigma + r(j)

         end do


*-----------------------------------------------------------------------
c for debug cKN
*-----------------------------------------------------------------------

      goto 990

      if( sigma .gt. 0.0d0 ) then
      sigmaa=0.d0

      do j = 1, 6
        sigmaa=sigmaa+r(j)
      end do

         prob = (sigma-sigmaa)/sigma

      if( ia .gt. 50 .and.
     &    prob .gt. 0.05 ) 
     &    write(6,'(''50-> '',2i3,e13.5)') ia,iz, prob

      if( ia .gt. 40 .and. ia .le. 50 .and.
     &    prob .gt. 0.05 ) 
     &    write(6,'(''40-> '',2i3,e13.5)') ia,iz, prob

      if( ia .gt. 30 .and. ia .le. 40 .and.
     &    prob .gt. 0.07 ) 
     &    write(6,'(''30-> '',2i3,e13.5)') ia,iz, prob

      if( ia .gt. 20 .and. ia .le. 30 .and.
     &    prob .gt. 0.3 ) 
     &    write(6,'(''20-> '',2i3,e13.5)') ia,iz, prob


      end if

*-----------------------------------------------------------------------

  990 continue

      return
      end  subroutine


************************************************************************
*                                                                      *
      subroutine fis(a,z,u,fisinh)
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally in LAHET code**************
C  FIS
c      Pick post fission parameters such as mass, charge, kinetic energy
c      and excitation energy
C=====================================================================
C <variables>
C     a    :   mass of nucleus before fission                     (IN)
C     z    :   charge  of nucleus before fission                  (IN)
C     u    :   excitation energy of nucleus before fission        (IN)
C  fisinh  : true=fission occur, false=no fission                 (IN)
C   ef     :   fission barrier                                    (IN)
C  zfis    :  charge of the fission fragment                     (OUT)
C  afis    :  mass of the fission fragment                       (OUT)
C  ufis    :  excitation energy of the fission fragment          (OUT)
C   er     :  recoil energy of the fission fragment              (OUT)
C  betf    :  recoil direction of the fission fragment           (OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)

      real *8rani, ranj, rijk, rijk0, ranb, rans
      common /rdata/ rijk, rijk0, rani, ranj, ranb, rans, nrand

      common /fiss2/ ef
      common /fiss3/ zfis(2),afis(2),ufis(2),er(2),betf(2,3)
      common /mom/ erec,bet0(3),ekin,bet1(3)

      common /options/alev, rcal, ifis

      parameter (amu=931.494d0)
      parameter (pi=3.1415926535898d0)

      logical fisinh
      dimension evodba(4)
c     following data for assymetric vs symetric picking
      data evodba/18.8d0,18.1d0,18.1d0,18.5d0/
      data sigmaa, aamean /6.5d0,140.d0/
c     asym gauss width,1.0/level density parameter and high mass mean.
      data afact, zfact, enmass /8.071323d0,0.782354d0,939.56563d0/
c     mass diff from 1 amu for neutron,diff p & n masses and n mass

      if (.not.fisinh) then
c***********************************************************************
c
c  logic error fissed called with fisinh unset.
c
c***********************************************************************
       write (*,60) a,z,u,erec
       return
      endif
c***********************************************************************
c
c       pick the masses
c
c***********************************************************************
c      z=real(jz)
c      a=real(ja)
      jz=nint(z)
      ja=nint(a)
      nck2=0
   10 continue
      nck2=nck2+1
      if (nck2.gt.10) then
c***********************************************************************
c
c     fission failure
c
c***********************************************************************
        write (*,50) ja,jz,u,erec,zfis(1),afis(1),zfis(2),afis(2)
        fisinh=.false.
        return
      endif
      temp=z*z/a
      if (temp.le.35.d0) then
        if (jz.le.88) go to 20
      elseif (u.le.62.d0) then
c***********************************************************************
c
c   high z fission mass distribution. competition for sym vs assym
c  simple symmetric to assymetric data fit
c
c***********************************************************************
        arg=-0.36d0*u
        arg=4.87d+03*exp(arg)
        proba=arg/(1.d0+arg)
c        if (fltrnf(dumm).le.proba) then
        if (rn(0).le.proba) then
c***********************************************************************
c
c      assymetric fission
c
c***********************************************************************
          a1=gaussn(aamean,sigmaa)
          go to 30
        endif
      endif
c***********************************************************************
c
c     find assymetric barrier for width computation
c     assymetric barrier from seaborg and vandenbosch
c     phys.rev. 88,507 (1952) phys.rev. 110,507 (1958)
c     find if ee eo oe or oo nucleus
c     na=1 odd-odd,2 even-odd,3 odd-even,4 even-even
c
c***********************************************************************
      if(a.gt.240) then        ! furi May 4, 1999
       in=ja-jz
       na=1
       if (jz.eq.2*(jz/2)) na=na+1
       if (in.eq.2*(in/2)) na=na+2
       temp=z*z/a
       ef=evodba(na)-0.36d0*temp
      endif                    ! furi May 4, 1999
   20 continue
C  furi ... May 4, 1999
      if(ifis.ne.0) then
       upr=u-ef
       if (upr.gt.1.d+02) upr=1.d+02
       sigmas=0.425d0*(upr*(1.d0-0.005d0*upr)+9.35d0)

      else

       xx=z*z/a
       bf=efms(z,a)
       if(bf.lt.0)bf=ef
cKN
       upr=min(u-bf,450.0d0)
c      upr = u - bf

       sigmas = 0.122d0*xx**2.d0 - 7.77 * xx + 134.0d0 + 3.32d-2 * upr

      endif
c***********************************************************************
c
c     sigmas is symmetrin fission mass width.taken from systematic of
c     neuzil & fairhall phys.rev. 129,2705,(1963)
c
c      low z fission is always symmetric
c
c      high z fission is sometimes.loop back here from 1 loop if
c      symmetric fission predicted.
c
c***********************************************************************
C  furi ... May 4, 1999
      amean=0.5d0*a
      a1=gaussn(amean,sigmas)
   30 continue
c***********************************************************************
c
c      1 loop for assymmetrin fission returns to here.
c
c***********************************************************************
      afis(1)=aint(a1)
c***********************************************************************
c
c    check for low final a
c
c***********************************************************************
      if (afis(1).lt.5.d0) afis(1)=5.d0
      if (a-afis(1).lt.5.d0) afis(1)=a-5.d0
      afis(2)=a-afis(1)
c***********************************************************************
c
c        pick the charge
c
c***********************************************************************
      z1=65.5d0*afis(1)/(131.d0+afis(1)**0.666667d0)
      z2=65.5d0*afis(2)/(131.d0+afis(2)**0.666667d0)
      z1=z1+.5d0*(z-z1-z2)
c***********************************************************************
c
c       we use constant charge density with a 2 unit gaussian smearing
c
c***********************************************************************
c   sigma = 0.75;   z1=gaussn(z1,dph)
c                                by furi 18/DEC/1997
c***********************************************************************
      if(ifis.ne.0) then
       sigz=2.0d0
      else
       sigz=0.75d0
      endif
      z1=gaussn(z1,sigz)
      zfis(1)=aint(z1)
      zfis(2)=z-zfis(1)
c***********************************************************************
c
c   check for reasonable z a combinations.
c
c***********************************************************************
      if (zfis(1).ge.afis(1)) go to 10
      if (zfis(2).ge.afis(2)) go to 10
      if (zfis(1).lt.1.d0) go to 10
      if (zfis(2).lt.1.d0) go to 10
c***********************************************************************
c
c      compute binding energy and actual masses of fragments
c
c***********************************************************************
      be0=afact*a-zfact*z-energm(jz,ja)
      rm0=enmass*a-zfact*z-be0
      iaf=nint(afis(1))
      izf=nint(zfis(1))
      be1=afact*afis(1)-zfact*zfis(1)-energm(izf,iaf)
      rm1=enmass*afis(1)-zfact*zfis(1)-be1
      iaf=nint(afis(2))
      izf=nint(zfis(2))
      be2=afact*afis(2)-zfact*zfis(2)-energm(izf,iaf)
      rm2=enmass*afis(2)-zfact*zfis(2)-be2
c***********************************************************************
c
c      pick recoil kinetic energy.use systematic of ...........
c      unik et.al. proc.3rd iaea symp.on phy.& chem. fision,rochester.vo
c
c***********************************************************************
      if(ifis.ne.0)then
       totkm=0.13323d0*z*z/a**0.33333333d0-11.4d0
      else
C... 4 May 1999 by furi
       x=z*z/a**0.33333d0
       if(x.le.900.d0) then
        totkm=0.131d0*x
       else if(x.le.1800.d0) then
        totkm=0.104d0*x+24.3d0
       else
        write(*,*)'Error in totk in subroutine fiss'
       endif
      endif
c***********************************************************************
c
c      use a width of 15% value at half height.
c
c***********************************************************************
      if(ifis.ne.0)then
       sigmak=0.084d0*totkm
C... 4 May 1999 by furi
      else
       if(x.lt.1000)then
        sigmak=86.5d0
       else if (x.lt.1800.d0)then
        sigmak=5.70d-4*(x-1000.d0)**2.d0+86.5d0
       else
        write(*,*)'Error in sigmak in subroutine fiss'
       endif
      sigmak=sqrt(sigmak)
      endif
c***********************************************************************
c
c   check event is energetically possible
c
c***********************************************************************
      temp2=u+be1+be2-be0
      nck=0
   40 continue
      totke=gaussn(totkm,sigmak)
      if (nck.gt.10) go to 10
      nck=nck+1
      if (totke.gt.temp2) go to 40
c***********************************************************************
c
c      pick excitation from equidistribution of original plus energy bal
c
c***********************************************************************
      temp=(temp2-totke)/a
      ufis(1)=afis(1)*temp
      ufis(2)=afis(2)*temp
cccc
c     find total masses, including excitation energies, at evap time
cccc
      amcf=rm0+u
      amc1=rm1+ufis(1)
      amc2=rm2+ufis(2)
      amdiff=amcf-amc1-amc2
      if (amdiff.lt.0.d0) go to 40
c***********************************************************************
c
c     amdiff= ekin should be satisfied
c
c***********************************************************************
C
c      ernff=erec
c      call fisdis

C... Velocity of pre-fission nucleus : vres*bet0
      pres=sqrt(erec*2+2.d0*a*amu*erec)
      vres=pres/(Erec+a*amu)

C... Momentum of CM system
      redm=amu*(afis(1)*afis(2))/a
      pcm=sqrt(totke**2+2.d0*redm*totke)
c      th  = acos( 2.0d0 * fltrnf(dumm) - 1.0d0 )
c      ph  = 2.0d0 * pi  * fltrnf(dumm)
      th  = acos( 2.0d0 * rn(0) - 1.0d0 )
      ph  = 2.0d0 * pi  * rn(0)

      pcmx = pcm* dsin(th) * dcos(ph)
      pcmy = pcm* dsin(th) * dsin(ph)
      pcmz = pcm* dcos(th)

C... Velocity of fission fragment 1 in CM system
      v1x=pcmx/sqrt(pcm**2 + (amu*afis(1))**2)
      v1y=pcmy/sqrt(pcm**2 + (amu*afis(1))**2)
      v1z=pcmz/sqrt(pcm**2 + (amu*afis(1))**2)

C... Velocity of fission fragment 2 in CM system
      v2x=-pcmx/sqrt(pcm**2 + (amu*afis(2))**2)
      v2y=-pcmy/sqrt(pcm**2 + (amu*afis(2))**2)
      v2z=-pcmz/sqrt(pcm**2 + (amu*afis(2))**2)

C... Boost to Lab system
      v1x=v1x+vres*bet0(1)
      v1y=v1y+vres*bet0(2)
      v1z=v1z+vres*bet0(3)

      v2x=v2x+vres*bet0(1)
      v2y=v2y+vres*bet0(2)
      v2z=v2z+vres*bet0(3)

C... Kinetic Energy and velocity of fission fragment 1 in Lab system
      gam=1.d0/sqrt(1.d0-v1x**2-v1y**2-v1z**2)
      er(1)= afis(1)*amu *(gam -1.d0)

      betf(1,1)=v1x/sqrt(v1x**2+v1y**2+v1z**2)
      betf(1,2)=v1y/sqrt(v1x**2+v1y**2+v1z**2)
      betf(1,3)=v1z/sqrt(v1x**2+v1y**2+v1z**2)

C... Kinetic Energy and velocity of fission fragment 2 in Lab system
      gam=1.d0/sqrt(1.d0-v2x**2-v2y**2-v2z**2)
      er(2)= afis(2)*amu *(gam -1.d0)

      betf(2,1)=v2x/sqrt(v2x**2+v2y**2+v2z**2)
      betf(2,2)=v2y/sqrt(v2x**2+v2y**2+v2z**2)
      betf(2,3)=v2z/sqrt(v2x**2+v2y**2+v2z**2)

      return
c
   50 format ('---> fission failed: ja=',i5,'  jz=',i5/'-       u=',
     1 1pe10.3,'      erec=',e10.3/'    zfis1=',e10.3,'     afis1=',
     2 e10.3/'     zfis2=',e10.3,'     afis2=',e10.3)
   60 format (//'  logic error in fiss.called with fisinh flag',' unse
     1t.',2i10,2f10.5)
      end subroutine


************************************************************************
*                                                                      *
C--------------------------------------------------------------------
      subroutine selectE(j,a,z,u,pekin,ekin)
C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      common /ejectl/ omega(70),ifa(70),ifz(70)
      common /emitr/ gj(70),q(70),V(70),delta(70),smalla(70)
      common /std1/ r(70),s(70),sigma,rr(70)
      common /fiss/ beta,alp

      common /sim/nocas

      aa=a-dble(ifa(j))
      izz=nint(z)-ifz(j)
      nn=nint(aa)-izz
      iaa=nint(aa)

      if(iaa.eq.1) then
       ekin=u-q(j)
       return
      endif

      bet=beta
      if(j.ne.1) bet = -V(j)

      ux = 2.5d0 + 150.d0 / aa
      ex = ux + delta(j)
      ax = getag(ex,izz,nn,isdum)
      tau = sqrt(ax / ux) - 1.5d0 / ux
      tau = 1.d0 / tau

      sx = 2.d0 * sqrt( ax * ux )
      e0 = ex - tau * (dlog(tau) - .25* dlog(ax) -1.25*dlog(ux) + sx)

      ppt=0.d0

      imax=1000

c.....Calculate excat decay width
c      tempt=0.0
c      tempt2=0.0
c      do 100 ifuri=0,imax
c        x=V(j)+(u-q(j)-V(j))/imax*ifuri
c        aaa=getag(u-q(j)-x,izz,nn,isdum)
c        temp=pe(u,q(j),delta(j),V(j),bet,tau,e0,ex,smalla(j),x)
c        temp2=pe(u,q(j),delta(j),V(j),bet,tau,e0,ex,aaa,x)
c
c        if(ifuri.eq.0)then
c         tempo=temp
c         tempo2=temp2
c         goto 100
c        endif
c
c        tempt=tempt+(tempo+temp)*(u-q(j)-V(j))/imax/2.d0*gj(j)
c        ddt=(tempo2+temp2)*(u-q(j)-V(j))/imax/2.d0*gj(j)
c        tempt2=tempt2+ddt
c        if(ddt/tempt2.lt.1.e-5) goto 101
c        tempo=temp
c        tempo2=temp2
c
c        tempm=tempt
c        write(44,*)x,tempt,tempt2
c 100  continue
c 101  continue
c      ran=pekin/r(j)/rr(j)
c      pekin=tempt2*ran

      do 30 ii=1,2
        do 10 i=0,imax
          x=V(j)+(u-q(j)-V(j))/imax*i

c          aaa=getag(u-q(j)-x,izz,nn,isdum)

          pp=pe(u,q(j),delta(j),V(j),bet,tau,e0,ex,smalla(j),x)

c          pp=pe(u,q(j),delta(j),V(j),bet,tau,e0,ex,aaa,x)

          if(i.eq.0)then
           ppo=pp
           goto 10
          endif

          dppt=(ppo+pp)*(u-q(j)-V(j))/imax/2.d0*gj(j)
          ppt=ppt+dppt
          ppo=pp
          ppm=ppt
c          write(55,*)i,x,ppm,pekin,j,r(j)*rr(j)

        if(ppm.gt.pekin) goto 20

 10   continue
      ran=pekin/r(j)/rr(j)
c      if(ran.gt..8d0) then
      if(ii.eq.1) then
       pekin=ppm*ran
       goto 30
      else

c      write(*,*)'error in selectE: Can not find Ekin at nocas= ',nocas
c      write(*,*)' jemiss=',j,', A=',a,', Z= ',z,', Ex/A= ',u/a,
c    &      ', random number =',ran
c      write(*,*)'pekin,ppm,rr,r', pekin,ppm,rr(j),r(j)
c      call parastop( 999 )
cKN
        ekin = -1.0
        return
cKN
      endif

C.....dubeg write
c       ppo=0.d0
c        do 11 iii=0,imax
c          x=V(j)+(u-q(j)-delta(j)-V(j))/imax*iii
c          pp=pe(u,q(j),delta(j),V(j),bet,tau,e0,ex,smalla(j),x)
c          if(iii.eq.0)then
c           ppo=pp
c           goto 11
c          endif
c          dppt=(ppo+pp)*(u-q(j)-delta(j)-V(j))/imax/2.d0*gj(j)
C          if(dppt/ppt.lt.1e-5) goto 40
c          if(dppt/ppt.lt.1e-5) write(*,*)'goto 40'
c          ppt=ppt+dppt
c          ppo=pp
c          ppm=ppt*exp(-smax)
c          write(*,*)iii,x,ppm,pekin,j,r(j)*rr(j)
c 11   continue

 30   continue

 20   ekin=x
c      write(70,*)x,pekin/r(j)/rr(j),ppm/tempm,j,i,tau
      return
      end subroutine

************************************************************************
*                                                                      *
C--------------------------------------------------------------------
      function pe(e,q,delta,V,bet,tau,e0,ex,smalla,x)
C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      if(e-q.le.0) then
       pe=0
       return
      endif

      if(x.gt.e-q-ex.and.x.le.e-q)then
c       write(*,*)'rho1',x,e-q-ex,e-q

cKN
       eag = (e-q-x-e0)/tau
c      eag = min( 200.d0, eag )
cKN
       pe = exp( eag ) / tau

      else if(x.ge.V.and.x.le.e-q-ex) then

cKN
       eag = 2.d0*sqrt(smalla*(e-q-delta-x))
c      eag = min( 200.d0, eag )
cKN
       pe = exp( eag )

       pe= pe/smalla**.25d0
       pe= pe/(e-q-delta-x)**1.25d0
      else
       pe=0.d0
      endif

c      write(*,*)pe,x,bet
      pe=pe*(x+bet)
      return
      end function

************************************************************************
*                                                                      *
C--------------------------------------------------------------------
      subroutine eye(j,aa,izz,nn,u,q,V,delta,smalla,beta,r,s)
C
C     variable       IN/OUT
C     j              I       emittor identifier
C     aa             I       residual mass after emission (=A-Aj)
C     izz            I       charge number of res nuclei after emission
C     nn             I       neutron number of res nuclei after emission
C     u              I       excitation energy of nuclei before emission
C     q              I       Q-value
C     V              I       kV in the equation
C     delta          I       pairing energy
C     smalla         I       level density parameter at U-Q-delta-kV
C     bet            I       beta in the equation
C     r              O       r in the equation
C                                      _______________
C     s              O       s = 2 \/a(U-Q-delta-kV)

C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

*-----------------------------------------------------------------------

            if( u - q - V .le. 0.0 .or. aa .le. 0.0 ) then

               r = 0.d0
               s = 0.d0
               return

            end if

*-----------------------------------------------------------------------

               ux  = 2.5d0 + 150.d0 / aa
               ex  = ux + delta

               ax  = getag(ex,izz,nn,isdum)

               tau = sqrt( ax / ux ) - 1.5d0 / ux
               tau = 1.d0 / tau
cKN
               t = ( u - q  - V ) / tau
               sx = 2.d0 * sqrt( ax * ux )

               ssx = sx - ex / tau
               ssx = min( 200.d0, ssx )

               eest =  exp( ssx ) * tau / ax**.25d0 / ux**1.25d0
cKN

         if( u - q - V .le. ex ) then

               eept = exp( t )

               eye1 = ( eept - 1.0 - t ) * tau
               eye0 =   eept - 1.0

               r = ( eye1 + ( beta + V ) * eye0 ) * eest

         else

cKN
               s  = 2.d0 * sqrt( smalla * ( u - q - delta - V ) )
               tx = ex / tau
cKN
               eeps = exp( s )
               eepx = exp( tx )
               sexp = exp( sx - s )

               eye0 = ( eepx - 1.0 ) * eest
               eye1 = (  ( eepx * ( t  - tx + 1.0 )
     &                 - ( t + 1.0 ) ) * tau ) * eest

               eye2 = ey2(s,sx,sexp)        * eeps
               eye3 = ey3(s,sx,smalla,sexp) * eeps


               r =  eye3 + eye1  + ( beta + V ) * ( eye0 + eye2 )

         end if

      return
      end subroutine


C--------------------------------------------------------------------
      function ey2(s,sx,sexp)
C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      func ( x ) =  1/x**1.5d0 + 1.5d0/x**2.5d0 + 3.75d0/x**3.5d0
c     &     + 13.125d0 * x**-4.5d0

      temp  = func(s)
      tempx = func(sx)

      ey2 = 2.d0 * sqrt(2.d0) *( temp - sexp * tempx )

      return
      end function

C--------------------------------------------------------------------
      function ey3(s,sx,a,sexp)
C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      ssqr = s**2
      sxqr = sx**2

      temp =     325.125d0 / s**4.5d0
     &     - sexp * ( ( 324.8d0   * ssqr + 3.28d0 * sxqr ) / sx**6.5d0 )
      temp = temp + 60.0d0 / s**3.5d0
     &     - sexp * ( ( 59.0625d0 * ssqr + 0.9375 * sxqr ) / sx**5.5d0 )
      temp = temp + 13.5d0 / s**2.5d0
     &     - sexp * ( ( 12.875d0  * ssqr + 0.625  * sxqr ) / sx**4.5d0 )
      temp = temp +  4.d0 / s**1.5d0
     &     - sexp * ( (  3.75d0   * ssqr + 0.25   * sxqr ) / sx**3.5d0 )
      temp = temp +  2.d0 / s**.5d0
     &     - sexp * ( (   1.5d0   * ssqr + 0.5    * sxqr ) / sx**2.5d0 )
      temp = temp
     &     - sexp * ( ( ssqr - sxqr ) / sx**1.5d0 )

      ey3 = temp / sqrt(2.d0) / a

      return
      end function

C--------------------------------------------------------------------
      function ey0(t)
C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      ey0 = exp(t)- 1.d0

      return
      end function

C--------------------------------------------------------------------
      function ey1(t,tx,tau)
C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      ey1 = ( exp(tx) * ( t  - tx + 1.d0 ) - ( t + 1.d0 ) ) * tau

      return
      end function

************************************************************************
*                                                                      *
      function energm(iz,ia)
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally in the LAHET code********
C  ENERGY
c    Calculate excess mass
C=====================================================================
C <variables>
C     ia   :   mass of nucleus                      (IN)
C     iz   :   charge  of nucleus                   (IN)
C  energm  :   Excess mass [MeV]                   (OUT)
C/////////////////////////////////////////////////////////////////////
      implicit doubleprecision(a-h,o-z)
      common /exims/ wapsm(0:150,0:250)

      logical isz, isn
      parameter (inn=150, izz=98)
      common /cook/ sz(izz), sn(inn), con(2), amean(240), pz(izz),
     1 pn(inn), isz(izz), isn(inn)

      data afact, zfact, enmass /8.071323d0,0.782354d0,939.56563d0/

      cam(a,z,b)=8.071323d0*a-0.782354d0*z-17.0354d0*a*(1.d0-1.84619d0
     1 *(a-2.d0*z)**2/a**2)+25.8357d0*b**2*(1.d0-1.71219d0*(a-2.d0*z)
     2 **2/a**2)*(1.d0-0.62025d0/b**2)**2+0.779d0*z*(z-1.d0)*(1.d0-1.5
     3 849d0/b**2+1.2273d0/a+1.5772d0/(a*b))/b-0.4323d0*exp(1.33333333
     4 33d0*log(z))*(1.d0-0.57811d0/b-0.14518d0/b**2+0.49597d0/a)/b

      energm=1.d10

      if(iz.eq.0.and.ia.gt.10) then
         energm=15.d0*ia
         return
      end if

      in=ia-iz
      if (iz.ge.0.and.iz.le.150.and.ia-iz.ge.0.and.ia-iz.le.250) then
         energm=wapsm(iz,ia-iz)
         if(energm.eq.0.d0) goto 10
         return
      endif

Cameron's mass formula
 10   if(iz.eq.6.and.ia.eq.12) return

      if( iz .eq. 0 ) return

      n=ia-iz
      a=dble(ia)
      z=dble(iz)
      a3=exp(0.3333333333d0*log(a))

cKN
      izk = min(iz,izz)
      ink = min(n,inn)
cKN
c////////////////&&&&&&&&&&
!   ink=           0  energm(iz,ia)
!   iz,ia =          10          10  n, inn=      0         150
!   iz,ia =          11          11  n, inn=      0         150
!    rather frequent for Fe 1500 MeV (KE/n) + p case
!
      if(ink <= 0 ) then
         write(0,*) ' ink=',ink, 'invalid in  energm(iz,ia)'
         write(0,*) 'iz,ia =', iz, ia, ' n, inn=',n, inn
         write(0,*) ' ink = 1 is used'
         ink = 1
      endif
c////////////////
      energm=cam(a,z,a3)+sz(izk)+sn(ink)

      return
      end function


************************************************************************
*                                                                      *
      function fprob (z,a,e)
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally in LAHET code**************
C  FPROB
C    Calculate fission probability
C=====================================================================
C <variables>
C     a    :   mass of nucleus before fission                     (IN)
C     z    :   charge  of nucleus before fission                  (IN)
C     e    :   excitation energy of nucleus before fission        (IN)
C  sigma   :   total decay width calculated in sub. gamma         (IN)
C   an     :   level density parameter for neutron emission
C                                calculated in sub. gamma         (IN)
C  frob    :   fission probability                               (OUT)
C   ef     :   fission barrier                                   (OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)
      parameter (inn=150, iiz=98)
      logical isz, isn
      common /cook/ sz(iiz), sn(inn), con(2), amean(240), pz(iiz),
     1 pn(inn), isz(iiz), isn(inn)

      common /std1/ r(70),s(70),sigma,rr(70)
      common /emitr/ gj(70),q(70),V(70),delta(70),smalla(70)
      common /fiss/ beta,alp
      common /fiss2/ ef

      real*8 i0,i1,hbarc
      parameter(amu=931.494d0,hbarc=197.327053d0)
c
c   function to compute the fission probability.
c     for z<90 uses ..........
c    uses statistical model fits.
      parameter (lnfis=18)
      dimension slope(lnfis), anort(lnfis)
      data slope /0.23d0,0.233d0,0.12225d0,0.14727d0,0.13559d0,0.15735d0
     1 ,0.16597d0,0.17589d0,0.18018d0,0.19568d0,0.16313d0,0.17123d0,
     2 6*0.17d0/
      data anort /219.4d0,226.9d0,229.75d0,234.04d0,238.88d0,241.34d0,
     1 243.04d0,245.52d0,246.84d0,250.18d0,254.0d0,257.8d0,261.3d0,
     2 264.8d0,268.3d0,271.8d0,275.3d0,278.8d0/

      data a1, a2, a3 /0.2185024d0,16.70314d0,321.175d0/
      data const /0.3518099d0/
      data c1, c2, c3 /1.089257d0,0.01097896d0,31.08551d0/

C   Initialization
      fprob=0.d0
      iz=nint(z)
      ia=nint(a)
      in=nint(a-z)
      u=e

      if (iz.gt.88) then
c
c     high z fission probability
c     for z> 89 and <100 ........
c     use the systematics of vandenbosch & huizenga with a ball park
c     observation that fission probability drops for most nuclei at 6 m

       iz=iz-88
       if (iz.gt.lnfis.or.u.lt.6.0) go to 20
       gamnf=slope(iz)*(a-anort(iz))
       gamnf=1.d1**gamnf
       agoes=1.d0
       go to 10
      endif

cKN
      izk = min(iz,iiz)
      ink = min(in,inn)

Calculate -Qn + paring energy + shell correction
c////////////////
      if(ink <= 1 ) then
         write(0,*) ' ink=',ink, 'fprob (z,a,e) '
         write(0,*) 'z,a,e =', z, a, e
         write(0,*) ' in, inn=',in, inn
         write(0,*) ' ink = 2 is used'
         ink = 2
      endif
c////////////////

      se=energm(iz,ia-1)+energm(0,1)-energm(iz,ia)-sz(izk)-sn(ink-1)
      if ((2*(iz/2).ne.iz).or.(2*(in/2).ne.in)) se=se-pz(izk)-pn(ink-1)
      if ((2*(iz/2).ne.iz).and.(2*(in/2).ne.in)) se=se-pz(izk)-pn(ink-1)
Calculate fission barrier
      x=z*z/a
      ef=x*(a1*x-a2)+a3+se
      if (ef.gt.u) go to 20    !Excited energy is below fission barrier
Calculate level density paameter for neutron emission
      an=0.125d0*(a-1.d0)
      an1=0.5d0/an
      an2=0.25d0*an1/an
Calculate level density parameter for fission
      af=x-c3
      af=an*(c1+c2*af*af)
      a1thrd=a**0.33333333d0
      ss=2.d0*sqrt(an*(u-se))
Calculate Gamma_n/Gamma_f
cKN
      if (ss.gt.10.d0) then
c     if (ss.gt.20.d0) then
cKN
C I0=J0
        i0=an1*(ss-1.d0)
C I1=J1
        i1=ss*an2*((ss+ss-6.d0)+6.d0)
        gamnf=const*(((0.76d0*i1-5.d-02*i0)*a1thrd+1.93d0*i1)*a1thrd+
     1       1.66d0*i0)
        s2=2.d0*sqrt(af*(u-ef))
        exps=0.d0
cKN
        eag = ss-s2
c       eag = min( 200.d0, eag )
cKN
        if (eag.gt.-150.d0) exps=exp(eag)

       gamnf=gamnf*exps*af/(s2-1.d0)
      else
        e1=exp(ss)
        i0=((ss-1.d0)*e1+1.d0)*an1
        i1=((6.d0+ss*(ss+ss-6.d0))*e1+ss*ss-6.d0)*an2
        gamnf=const*(((0.76d0*i1-5.d-02*i0)*a1thrd+1.93d0*i1)*a1thrd+
     1  1.66d0*i0)
        ss=2.d0*sqrt(af*(u-ef))
        e2=((ss-1.d0)*exp(ss)+1.d0)/af
        gamnf=gamnf/e2
      endif
C....Since we could not find any reference of the following suppression
C     factor for high excited energy fission, i.e. agoes, we did not
C     mention agoes in the paper.
C     But we include agoes because it is used in LAHET code.
      call drein1 (1,s(1),smalla(1),eye1,eye0)
      call drein2 (s(1),smalla(1),eye2)
      if( eye1*eye0 .ne. 0.0d0 ) then
      epsav=(eye2+beta*eye1)/(eye1+beta*eye0)
      agoes=(u-7.d0)/(epsav+7.d0)
      else
      agoes = 1.0
      end if
      if (agoes.lt.1.d0) agoes=1.d0
      goto 10

 10   continue
      fprob=1.d0/(1.d0+gamnf)/agoes
 20   continue
      return
      end function

************************************************************************
*                                                                      *
      function dostg(i,z)

C=====================================================================
C     This routine is originally in HETC code
C=====================================================================

      implicit doubleprecision(a-h,o-z)

      dimension t(3,4)

*-----------------------------------------------------------------------
*     Set Dostrovsky's parameter set,
*     the footnote of PR116(1959)683 (p.699)
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------

C      i=1 -> Calculate kp
C      i=2 -> Calculate k_alpha
C      i=3 -> Calculate cp

*-----------------------------------------------------------------------
*        kp
*-----------------------------------------------------------------------

      data ( t(1,i), i = 1, 4 ) /
     &      0.51d0,    !  t(1,1)
     &      0.60d0,    !  t(1,2)
     &      0.66d0,    !  t(1,3)
     &      0.68d0/    !  t(1,4)

*-----------------------------------------------------------------------
*        k_alpha
*-----------------------------------------------------------------------

      data ( t(2,i), i = 1, 4 ) /
     &      0.81d0,    !  t(2,1)
     &      0.85d0,    !  t(2,2)
     &      0.89d0,    !  t(2,3)
     &      0.93d0/    !  t(2,4)

*-----------------------------------------------------------------------
*        cp
*-----------------------------------------------------------------------

      data ( t(3,i), i = 1, 4 ) /
     &       0.0d0 ,   !  t(3,1)
     &      -0.06d0,   !  t(3,2)
     &      -0.10d0,   !  t(3,3)
     &      -0.10d0/   !  t(3,4)

*-----------------------------------------------------------------------

      if (z-50.d0) 30,10,10
   10 dostg=t(i,4)
   20 return
   30 if (z-20.d0) 40,40,50
   40 dostg=t(i,1)
      go to 20
   50 n=.1d0*z
      x=10.d0*(n+1.d0)
      x=(x-z)*.1d0
      dostg=x*t(i,n-1)+(1.d0-x)*t(i,n)
      go to 20
      end function


************************************************************************
*                                                                      *
C********This routine is the same as the drein1 in LAHET code**********
      subroutine drein1 (j,s,a,eye1,eye0)
      implicit doubleprecision(a-h,o-z)
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,
     2 dp2th=dp2/dp3)
c
c     compute statistical theory emission integrals
c     for s<dph use a series expansion.
c     for s>dph the explicit relationship
c     return for neutrons (and compute eye0)
c     return 1 for all others.
c
c     coeficients for series expansions
c
c     correction of c1 by r. e. prael
      dimension c0(7), c1(7)
      data c0 /0.66666667d0,0.25d0,0.06666667d0,0.01388889d0,
     1 0.00238095d0,0.00034722d0,0.00004409d0/
      data c1 /0.53333333d0,0.16666667d0,0.03809524d0,0.00694444d0,
     1 0.0010582d0,0.00013889d0,0.0000160d0/
c
      exps=dp0
      if (s.lt.1.d+02) exps=exp(-s)
      if (s.lt.dph) go to 10
c///// explicit relation
      b=dph/a
      eye1=b*b*(dp3+s*(s-dp3)+exps*(dph*s*s-dp3))
      if (j.eq.1) eye0=b*(s-dp1+exps)
      return
c///// small s series expansion
   10 continue
c
c  eye1=(1.0/(8*a*a))*(s**4/4)*(sum n=0 to 7:8*s**n/(n!*(n+2)*(n+4))
c
      eye1=dp1
      b=dp1
      do 20 n=1,7
      b=b*s
      c=b*c1(n)
      if (c.lt.1.0d-7) go to 30
      eye1=eye1+c
   20 continue
   30 continue
      b=s*s/a
      eye1=eye1*exps*b*b*0.03125d0
      if (j.gt.1) return
c eye0 (neutrons only)=(.5/a)*s**2/2*(sum n=0 to 7:2*s**n/(n!*(n+2))
      eye0=dp1
      b=dp1
      do 40 n=1,7
      b=b*s
      c=b*c0(n)
      if (c.lt.1.0d-7) go to 50
      eye0=eye0+c
   40 continue
   50 continue
      eye0=exps*eye0*s*s*0.25d0/a
      return
      end subroutine


************************************************************************
*                                                                      *
C*********This routine is the same as the drein2 in LAHET code**********
      subroutine drein2 (s,a,eye2)
      implicit doubleprecision(a-h,o-z)
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,
     2 dp2th=dp2/dp3)
c
c     compute statistical theory emission integrals
c
c     compute third integral
c
c     for s<dph use a series expansion.
c     for s>dph the explicit relationship
c
c     coeficients for series expansions
c
      dimension c2(7)
      data c2 /0.45714286d0,0.125d0,0.02539683d0,0.00416667d0,
     1 0.0005772d0,0.00006944d0,0.0000074d0/
      exps=dp0
      if (s.lt.1.d+02) exps=exp(-s)
      if (s.lt.dph) go to 10
c///// explicit relation
      b=s*s
      eye2=0.25d0*(s*(15.d0-s*(6.d0-s))-15.d0+(15.d0+0.125d0*b*(b-12.d0)
     1 )*exps)/(a*a*a)
      return
   10 continue
c
c///// series expansion
c    eye2=(1/(32a**3))*s**6/6*(sum n=0 to 7:48*s**n/(n!(n+2)(n+4)(n+6))
c
      eye2=dp1
      b=dp1
      do 20 n=1,7
      b=b*s
      c=b*c2(n)
      if (c.lt.1.0d-7) go to 30
      eye2=eye2+c
   20 continue
   30 continue
      b=0.25d0*s*s/a
      eye2=eye2*b*b*b*exps*0.33333333d0
      return
      end subroutine 


************************************************************************
*                                                                      *
      function efms(z,a)
C/////////////////////////////////////////////////////////////////////
C  EFMS
C  Fission barrier given by Myer & Swaiteski (PRC60,014606,1999)
C=====================================================================
C <variables>
C     a   :   the mass of a fissioning nucleus      (IN)
C     z   :   the charge of a fissioning nucleus    (IN)
C   efms  :   fission barrier  [MeV]                (OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)

      common /fiss4/ shell(150,250)

C... 8/15/1999
      parameter(x0=48.5428d0, x1=34.15d0)
      f1(t)=1.99749d-4*(x0-t)**3
      f2(t)=5.95553d-1-0.124136*(t-x1)

      efms=0.d0
      iz=nint(z)
      ia=nint(a)

      c=1.9+(z-80)/75
      ai=1.-2*(z/a)
      xx=1-c*ai**2
      ss=a**.66667*xx
      x=z**2/a/xx
      if(ia-iz.le.0.or.ia-iz.gt.250.or.iz.gt.150.or.iz.lt.1) then
       sh=0.d0
      else
       sh=shell(iz,ia-iz)
      endif
      if(x.ge.x1.and.x.le.x0) then
       efms=ss*f1(x)-sh
      else if(x.ge.20.and.x.lt.x0)then
       efms=ss*f2(x)-sh
      else
       efms=-1.0
      endif

      return
      end function


************************************************************************
*                                                                      *
       function gaussn (xmean,sd)
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally from the LAHET code********
C  GAUSSN
c      Gaussian randum number gemerator
c  compute random gaussian number for given
c  mean and s.d.
c  uses mean of sum of 12 uniform r.n"s
c
C/////////////////////////////////////////////////////////////////////
      implicit doubleprecision(a-h,o-z)

      a=0.d0
      do 10 n=1,12
        a=a+rn(0)
 10   continue
      gaussn=(a-6.d0)*sd+xmean
      return
      end function

************************************************************************
*                                                                      *
      function getag(e,iz,in,is)
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally in LAHET code**************
C  GETA
c      Calculate level density parameter
C=====================================================================
C <variables>
C     iz   :   charge  of nucleus                                 (IN)
C     in   :   neutron number of nucleus                          (IN)
C     e    :   excitation energy of nucleus                       (IN)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)

      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,
     2 dp2th=dp2/dp3)

      logical isz, isn
      parameter (inn=150, izz=98)
      common /cook/ sz(izz), sn(inn), con(2), amean(240), pz(izz),
     1 pn(inn), isz(izz), isn(inn)
      common /options/alev, rcal, ifis

      ia=in+iz
      aa=ia

      if(alev.eq.0.) go to 10

C     a=A/a0 ...simple level density parameter
      if(alev.eq.1.d0) then
       getag=aa/8.d0
      else
       getag=aa/alev
      endif
      return

C Gilbert-Cameron-Cook-Ignatyuk(GCCI) level density parameter

 10   if (iz.gt.izz.or.in.gt.inn) then
cKN
         getag = aa / 8.0d0
         return

c      write(*,*)'ERROR: Can not calculate the GCCI level density
c    & parameter.  IZ and IN are ',iz,in
c      call parastop( 999 )

      endif

      if (in.ge.9.and.iz.ge.9) then
        if (isz(iz).or.isn(in)) then
          is=2
        else
          is=1
        endif
        st=sz(iz)+sn(in)
        fa=(9.17d-3*st+con(is))
      else
        if (ia.le.25) then
         fa=0.125d0
        else if (ia.gt.25.and.ia.le.240) then
         fa=amean(ia)/aa
        else
         fa=(amean(240)+(2.5d0-0.1d0*amean(240))*(real(ia)-240.d0))/aa
        endif
      endif

      deldel=paire(iz,in)
      u=max(5.d-02*(e-deldel),1.d-05)
      gu=(dp1-exp(-u))/u
      temp=aa*(fa*gu+(dp1-gu)*(0.1375d0-8.36d-05*aa))
cKN
      getag = abs( temp )

      return
      end function


************************************************************************
*                                                                      *
      function paire(iz,in)

      implicit doubleprecision(a-h,o-z)

      logical isz, isn
      parameter (inn=150, iiz=98)
      common /cook/ sz(iiz), sn(inn), con(2), amean(240), pz(iiz),
     1 pn(inn), isz(iiz), isn(inn)

cKN
      jz = min(iz,iiz)
      jn = min(in,inn)
cKN

      if(in.eq.0) jn=1
      if(iz.eq.0) jz=1
      paire=pz(jz)+pn(jn)
c      write(*,*)in,jn,iz,jz,paire
      return

      end function


************************************************************************
*                                                                      *
      function rb(a,ia,j)
C/////////////////////////////////////////////////////////////////////
C  RB
C    Calculate Nuclear radius for a geometric cross section
C=====================================================================
C <variables>
C     a   :   mass of nucleus #1                    (IN)
C     z   :   charge  of nucleus #1                 (IN)
C    ia   :   mass of nucleus #2                    (IN)
C    iz   :   charge  of nucleus #2                 (IN)
C     j   :   type of the nufleus #2                (IN)
C    ck   :   transmission probability              (IN)
C   voul  :   Coulomb potential  [MeV]              (OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)

      common /options/alev, rcal, ifis
      parameter (r0=1.5d0)

      if(rcal.gt.0.0) goto 20

C... Dostrovsky et. al.

      if (j.le.6) then


         if( j .gt. 2 ) then

            rb = a**.333333d0 + dble(ia)**.333333d0

         else

            rb = a**.333333d0

         end if

       rb = rb * r0

C... Matsuse et al. PRC26(1982)2338

      else

       r1=1.12d0*a**.333333d0-0.86d0/a**.333333d0
       r2=1.12d0*dble(ia)**.333333d0-0.86d0/dble(ia)**.333333d0
       rb=r1+r2+2.85d0

      endif

       return

C...Simple form

 20   continue

      if(rcal.eq.10.0)then
       rr=1.5d0
      else
       rr=rcal
      endif
      r1=rr*a**.333333d0
      r2=rr*dble(ia)**.333333d0
      rb=r1+r2
      return
      end function


************************************************************************
*                                                                      *
C--------------------------------------------------------------------
      function rho(ia,iz,u,delta,smalla)
C
C     variable       IN/OUT
C     ia             I       residual mass after emission (=A-Aj)
C     iz             I       residual charge after emission 
C     u              I       excitation energy of nuclei before emission
C     delta          I       pairing energy
C     smalla         I       level density parameter at U-Q-delta-kV

C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

         nn = ia - iz
         ux = 2.5d0 + 150.d0 / dble(ia)
         ex = ux + delta

      if( u .ge. ex ) then
cKN
         eag =  2.d0 * sqrt( smalla * ( u - delta ) )
c        eag =  min( 200.d0, eag )
cKN
         rho = exp( eag )
     &       / smalla**.25d0 / ( u - delta )**1.25d0

      else

         ax  = getag(ex,iz,nn,isdum)
         tau = 1.0 / ( sqrt( ax / ux ) - 1.5d0 / ux )
cKN
         eag = 2.d0 * sqrt( ax * ux ) - ( ex - u ) / tau
c        eag =  min( 200.d0, eag )
cKN
         rho = exp( eag )
     &       / ax**.25d0 / ux**1.25d0

      end if

       if(rho.le.0.d0) write(*,*)'error in GEM, rho =',rho

      return
      end function


************************************************************************
*                                                                      *
      function vcoul(z,a,iz,ia,ck,j)
C/////////////////////////////////////////////////////////////////////
C  VCOUL
C    Calculate Coulomb potential
C=====================================================================
C <variables>
C     a   :   mass of nucleus #1                    (IN)
C     z   :   charge  of nucleus #1                 (IN)
C    ia   :   mass of nucleus #2                    (IN)
C    iz   :   charge  of nucleus #2                 (IN)
C     j   :   type of the nufleus #2                (IN)
C    ck   :   transmission probability              (IN)
C   voul  :   Coulomb potential  [MeV]              (OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)

      parameter (rc=1.70d0, ee=137.0359895d0, hbarc=197.327053d0)
      common /options/alev, rcal, ifis


C  No coulomb potential for neutron emission

      if(j.eq.1) then
       vcoul=0.d0
       return
      endif


      if(rcal.gt.0.0) goto 30

      if(j.le.6) then

C Dostrovsky's parameter set

       vcoul=hbarc/rc/ee*ck
       r2 = dble(ia)**.333333d0

       if(j.le.6) r2=1.2d0/rc
       if(j.eq.2) r2=0.d0

       r1 =a**.333333d0
       if(a.le.4.and.z.le.2) r1=1.2d0/rc

       r0=r1+r2

       goto 20

      else

C Matsuse's parameter set ...PRC26(1982)2338

       vcoul=hbarc/ee

       r1=1.12d0*a**0.333333d0-0.86d0/a**0.333333d0
       r2=1.12d0*dble(ia)**.333333d0-0.86d0/dble(ia)**.333333d0

       r0=r1+r2+3.75d0

       goto 20

      endif

C...Simple parameter set
 30   vcoul=hbarc/ee
      if(rcal.eq.10.0) then
       rr=1.5d0
      else
       rr=rcal
      endif
      r1=rr*a**.333333d0
      r2=rr*dble(ia)**.333333d0
      r0=r1+r2

 20   vcoul = vcoul * dble(iz) * z / r0

      return
      end function
