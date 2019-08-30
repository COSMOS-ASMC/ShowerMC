c  f  ajmrr2   to calculate relativistic 2-body distance.
c  f  ajmpp2   to calculate relativistic 2-body distance.
c***********************************************************************

c...This is under construction!
      subroutine jamclust

c...Purpose: to determine nuclear cluster.
      include 'jam1.inc'
      include 'jam2.inc'

c...Maximum baryon number.
      parameter(mxb=500)
      dimension mascl(mxb),num(mxb),it(0:mxb),itc(0:mxb,0:mxb)
      dimension isort(mxb),isorti(mxb)
      dimension rhoa(mxb)
      dimension ntb(10)
c...ntb(1) : neutron number of the cluster
c...ntb(2) : proton number of the cluster
c...ntb(3) : lambda number of the cluster
c...ntb(4) : sigma-  number of the cluster
c...ntb(5) : sigma0  number of the cluster
c...ntb(6) : sigma+  number of the cluster
c...ntb(7) : xi-  number of the cluster
c...ntb(8) : xi0  number of the cluster
c...ntb(9) : Omega-  number of the cluster
c...ntb(10): Other(resonance) number of the cluster
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      mstd(91)=mstd(91)+1

c....Transport all particles.
      if(pard(1).lt.1000) then
        call trsporta(1000.0D0)
      endif

c...Defiene cluster distance.
      rcc2=parc(151)**2
      pcc2=parc(152)**2

c...Find boost parameters becase output is Lab.system.
c     if( mstc(4) .eq. 0 ) then
c       betac = 0.0
c       gamac = 1.0
c     else
c       betac = pard(5)
c       gamac = pard(6)
c     end if
 
c     cpf2 = ( 1.5 * paru(1)**2 * 
c    $    ( 4.0 * paru(1) * parc(15))**(-1.5) )**(2./3.)
c    &     * paru(3) ** 2
c
c...Calculate overlap of wave packet.
c     c0w=0.25/parc(15)
c     do i = 1, nbary
c       rhoa(i) = 0.0
c       do j = 1, nbary
c         rhoa(i) = rhoa(i) + rha(i,j)
c       end do
c       rhoa(i) = ( rhoa(i) + 1.0 ) ** (1./3.)
c     end do
c
c     do i=1,nbary
c       rhoa(i)=0.0
c     do j=1,i-1
c     do j=1,i
c       dns=exp(-ajmrr2(i,j)*c0w)
c       rhoa(i)=rhoa(i)+dns
c       rhoa(j)=rhoa(j)+dns
c     end do
c     end do

c...Identify the nuclear clusters.
      do i=1,mxb
        mascl(i)=1
        num(i)=i
      end do
      nclst = 1
      ichek = 1
      do i=1,nbary-1
        j1=ichek+1
        do j=j1,nbary
          rdist2=ajmrr2(num(i),num(j))
          pdist2=ajmpp2(num(i),num(j))
c         pcc2=cpf2*(rhoa(num(i))+rhoa(num(j)))**2
c....Cluster
          if(rdist2.lt.rcc2.and.pdist2.lt.pcc2) then
            ibuf=num(ichek+1)
            num(ichek+1)=num(j)
            num(j)=ibuf
            ichek=ichek+1
            mascl(nclst)=mascl(nclst)+1
          end if
        end do
        if(ichek.eq.i) then
          nclst=nclst+1
          ichek=ichek+1
        end if
      end do

c...Sort for summary.
      do i=1,nclst
        isort(i)=i
      end do

 510  continue
      nexch = 0
      do i=1, nclst-1
        if(mascl(isort(i)).lt.mascl(isort(i+1))) then
          nexch=nexch+1
          ibuf=isort(i)
          isort(i)=isort(i+1)
          isort(i+1)=ibuf
        end if
      end do
      if(nexch.ne.0) goto 510

      do i=1,nclst
        isorti(isort(i))=i
      end do

      itc(0,0) = nclst
      inum = 0
      do i=1,nclst
        itc(isorti(i),0)=mascl(i)
        do j=1,mascl(i)
          inum=inum+1
          itc(isorti(i),j)=num(inum)
        end do
      end do


      nclust=0  ! number of clusters
c=====================================================================*
      do 1000 i = 1, nclst       !...Loop over all cluster.
c=====================================================================*

        mclst=itc(i,0)

c...This is a baryon.
        if(mclst.eq.1) goto 1000

        it(0) = itc(i,0)
        ip=itc(i,1)

        if(mclst.eq.2) then
          kf1=k(2,itc(i,1))
          kf2=k(2,itc(i,2))
          if(kf1.eq.2212.and.kf2.eq.2212) goto 1000 ! pp bound state?
          if(kf1.eq.2112.and.kf2.eq.2112) goto 1000 ! nn bound state?
c....deuteon: spin factor.
          if((kf1.eq.2212.and.kf2.eq.2112)
     $     .or.(kf1.eq.2112.and.kf2.eq.2212)) then
             if(rn(0).gt.3.D0/4.D0) goto 1000
          endif
        endif

c....Momentum and kinetic energy of cluster.
        pclx=0.0D0
        pcly=0.0D0
        pclz=0.0D0
        esum=0.0D0
        rclx=0.0D0
        rcly=0.0D0
        rclz=0.0D0
        do l=1,10
          ntb(l)=0
        end do

c......Loop over mass of one cluster,determine property of cluster.
c-------------------------------------------------
        do j = 1, itc(i,0)
c-------------------------------------------------
            it(j)=itc(i,j)
            jp=it(j)
            k(1,jp)=11
            kf=k(2,jp)
            if( kf .eq. 2112) then 
              ntb(1)=ntb(1)+1
            else if(kf.eq.2212) then
              ntb(2)=ntb(2)+1
            else if( kf .eq. 3122) then 
              ntb(3)=ntb(3) + 1
            else if( kf .eq. 3112 ) then 
              ntb(4)=ntb(4) + 1
            else if( kf .eq. 3212 ) then 
              ntb(5)=ntb(5) + 1
            else if( kf .eq. 3222 ) then 
              ntb(6)=ntb(6) + 1
            else if( kf .eq. 3312) then 
              ntb(7)=ntb(7) + 1
            else if( kf .eq. 3322) then 
              ntb(8)=ntb(8) + 1
            else if( kf .eq. 3334 ) then 
              ntb(9)=ntb(9) + 1
            else
              ntb(10)=ntb(10) + 1
            end if

             pclx = pclx + p(1,jp)
             pcly = pcly + p(2,jp)
             pclz = pclz + p(3,jp)
             esum = esum + p(4,jp)
             rclx = rclx + r(1,jp)
             rcly = rcly + r(2,jp)
             rclz = rclz + r(3,jp)
c-------------------------------------------------
        end do  !*** end loop over one cluster
c--------------------------------------------------

c...Save cluster information to the arrays.
c....cmrest   : Total mass of the cluser per nucl..
c....excit    : Excitation energy of the cluster per nucl.  (MeV)
c....jj       : Angular momentum of the cluser.

        nclust=nclust+1
        jj=0
        excit=0.0D0
        p(1,ip)=pclx
        p(2,ip)=pcly
        p(3,ip)=pclz
        p(4,ip)=esum
        p(5,ip)=excit
        k(1,ip)=5
        k(2,ip)=1000000*ntb(1)+1000*ntb(2)+ntb(3)+1000000000
        k(3,ip)=1000000*ntb(4)+1000*ntb(5)+ntb(6)+1000000000
        k(4,ip)=1000000*ntb(7)+1000*ntb(8)+ntb(9)+1000000000
        k(5,ip)=ntb(10)
        k(6,ip)=jj
        k(7,ip)=0
        k(8,ip)=0
        k(9,ip)=mclst  ! baryon number
        k(10,ip)=0
        k(11,ip)=0
        r(1,ip)=rclx/mclst
        r(2,ip)=rcly/mclst
        r(3,ip)=rclz/mclst

c=======================================================================
 1000 continue  !............... end of cluster loop.
c=======================================================================

c...Save cluster number.
      mstd(92)=nclust

c...Remove unwanted entries.
      call jamedit
 
      end

c**************************************************** 

      function ajmrr2(j,i)

c...Purpose to calculate relativistic 2-body distance.
      include 'jam1.inc'
      include 'jam2.inc'

        eij  = p(4,j) + p(4,i)
        drx  =   r(1,j) - r(1,i)
        dry  =   r(2,j) - r(2,i)
        drz  =   r(3,j) - r(3,i)
        dpx  =   p(1,j) - p(1,i)
        dpy  =   p(2,j) - p(2,i)
        dpz  =   p(3,j) - p(3,i)

        dbx  = ( p(1,j) + p(1,i) ) /  eij
        dby  = ( p(2,j) + p(2,i) ) /  eij
        dbz  = ( p(3,j) + p(3,i) ) /  eij

        bij2 = dbx*dbx + dby*dby + dbz*dbz
        gam2  = eij**2
     $        /(eij**2
     $    -((p(1,j)+p(1,i))**2+(p(2,j)+p(2,i))**2+(p(3,j)+p(3,i))**2))

        ajmrr2  = drx*drx + dry*dry + drz*drz
     $          + gam2*(drx*dbx+dry*dby+drz*dbz)**2


      end

c**************************************************** 

      function ajmpp2(j,i)

c...Purpose to calculate relativistic 2-body momentum distance.
      include 'jam1.inc'
      include 'jam2.inc'

        dpx  =   p(1,j) - p(1,i)
        dpy  =   p(2,j) - p(2,i)
        dpz  =   p(3,j) - p(3,i)
        eij=p(4,i)+p(4,j)
        dbx  = ( p(1,j) + p(1,i) ) /  eij
        dby  = ( p(2,j) + p(2,i) ) /  eij
        dbz  = ( p(3,j) + p(3,i) ) /  eij
        bij2 = dbx*dbx + dby*dby + dbz*dbz
        gam2  = eij**2/(eij**2
     $    -((p(1,j)+p(1,i))**2+(p(2,j)+p(2,i))**2+(p(3,j)+p(3,i))**2))
        ajmpp2  = dpx*dpx + dpy*dpy + dpz*dpz
     $      +( -(p(4,j) - p(4,i))**2
     $      + gam2 * ( ( p(5,j)**2 - p(5,i)**2 ) / eij )**2 )

      end
