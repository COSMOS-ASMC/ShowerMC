!      mugeo0: init for coord. transformation
!      muptos: 1ry to *
!      muroc0: init. rock profile
!      murock: compute path length for given zenith and fai in rock
!      muradl: compute r.l of rock or other whose z, z'**2/a are given
!           rndc    !  not included
!           kmover  !
!           kcossn  !
!           kgauss  !subroutines to be imported from kklib
!           k4ptdi  !
!           kmover  !
!     ****************************************************************
!     *                                                              *
!     * mugeo0: initialize routine for computing geomagnetic field   *
!     *         effect (give gmf constants and 1ry direction)        *
!     *                                                              *
!     * muptos: 1ry to * system conversion                           *
!     *                                                              *
!     ************************ tested 90.02.23 ***********************
!
!  /usage/   call mugeo0(bh, bv, w01,  w02, w03, t)
!                 this must be called when the 1ry direction is fixed
! coordinate system:
!
!          let y* be geomagnetic north
!              x*                west
!              z* vertical axis going downwards
!
!      then  gmf is // to y*-z* plane
!
!          let the unit vector of 1ry paticle be z#=(w01, w02, w03) in
!          * system.  ( # is for vector)
!          form x axis by  b# x z#  where b# is the vector of gmf so
!          that  x# = b# x z# / ! b# x z# ! (if b#//z#, x should be
!          x* ).    y axis is on b#-z# plane, i.e y# = z# x x#.
!          the angle between b# and z# is alfa.   gmf effect is
!          computed easily in the frame where x axis is rotated
!          so that z axis coincide with b#.    then the deflection etc
!          is transformed to (x,y,z) system where coordinate of ptcls
!          are measured.
!
!
!  *** mugeo0 ***
!
! -- input --
!     bh:  horizontal component of gmf             (gauss)
!     bv:  vertical   //               (may be < 0)  //
!    w01:  1ry ptcle direction cosine to x* axis
!    w02:  1ry ptcle direction cosine to y* axis
!    w03:  1ry ptcle direction cosine to z* axis
!
! -- output --
!c     t:  transformaton matrix t(3,3).
!
!               x*    {t(1,1),t(1,2),t(1,3)}  x
!               y* =  {t(2,1),t(2,2),t(2,3)}  y
!               z*    {t(3,1),t(3,2),t(3,3)}  z
!
!
      subroutine mugeo0(bh, bv, w01, w02, w03, t)
!
      dimension t(3,3)
      data small /1.e-3/
      dimension tm(3,3)
!
!          mag. of gmf
      b=sqrt(bh**2 + bv**2)
!          cos(alfa)
      cosa=(bh * w02  +  bv * w03)/b
      cosa2=cosa**2
      sina=sqrt(1. - cosa2)
      bt=b * sina
      sb2= w01**2 + (bh/b-w02)**2 + (bv/b-w03)**2
      if(sb2 .gt. small) then
          cota=cosa/sina
          tm(1,1)= ( bh*w03 - bv*w02 ) /bt
          tm(1,2)= -cota*w01
          tm(1,3)= w01
          tm(2,1)= bv*w01/bt
          tm(2,2)= bh/bt - cota*w02
          tm(2,3)= w02
          tm(3,1)= -bh*w01/bt
          tm(3,2)= bv/bt - cota*w03
          tm(3,3)= w03
      else
          tm(1,1)= 1.
          tm(1,2)= 0.
          tm(1,3)= 0.
          tm(2,1)= 0.
          tm(2,2)=bv/b
          tm(2,3)=bh/b
          tm(3,1)=0.
          tm(3,2)= -bh/b
          tm(3,3)= bv/b
      endif
      call kmover(tm, 1, 9, t, 1)
      return
!
!       convert (x,y,z) in 1ry system into * system
!     ************
      entry muptos(x, y, z, xs, ys, zs)
!     ************
!
      tmp1= tm(1,1)*x + tm(1,2)*y + tm(1,3)*z
      tmp2= tm(2,1)*x + tm(2,2)*y + tm(2,3)*z
      zs  = tm(3,1)*x + tm(3,2)*y + tm(3,3)*z
      xs  = tmp1
      ys  = tmp2
      return
      end
!c        test muroc0; murock
!c
!        parameter (ipdb=99, pi=3.141592, Torad=pi/180.)
!        parameter (nr=61, ntet=360)
!        dimension dh(nr, ntet)
!        dimension tet(ntet), r(nr)
!        character*8  picnm
!c
!c       call muroc0('c2g5100.rock.data(kgf)', a, z, zba, z2ba,rho, bh,
!c   *   bv, beta)
!c       call muroc0('c2g5100.rock.data(frejus)', a, z, zba, z2ba,rho,
!c   *   bh, bv, beta)
!        dz=1.
!        da=1.
!        do 100 i=1, ntet
!            t=(i-1)*Torad
!            do 90 j=1, nr
!                 rt=(j-1)*Torad
!                 call murock(rt, t, dh(j, i) )
!  90        continue
! 100    continue
!c
!c         this is to draw  (length for given teta, fai)
!        call opnpdb('c2s5001.#pdb.data&', ipdb)
!                 do 30 i=1, ntet
!                    tet(i)=((i-1) )*Torad
!  30             continue
!                 do 40 i=1, nr
!                    r(i)=  (i-1)
!  40             continue
!                 picnm='frjtf'
!                 call newpic(picnm, 1,  'length&')
!                 call cont2c(nr, ntet, r, tet, dh, nr)
!                 call picend
!        call clspdb
!     end
!    ***********************************************************
!    *
!    * muroc0: init. for geometry of the detector position underground
!    * murock: comute path length of muon to reach the detector
!    *
!    ***********************************************************
!
!
!/usage/   call muroc0(dsn, a, z, zba, z2ba, rho, bh, bv, beta)
!          call murock(tet, fai, path)
!    dsn: input. dataset name which contatins path length table for
!                given teta and fai.
!                for the sturucter of the data, see
!                'c2g5100.rock.data(frejus)'.  or ..(kgf)
!     a,z,zba,z2ba,rho:  output. rho is in g/cm**2
!           bh, bv: output. geomagnetic strength in gauss. horizontal
!                           and vertical component.
!         beta: angle btween x* and detector direction (in deg).
!
!    tet: input. zenith angle (====in rad===)
!    fai: input. azimuthal angle (in rad).  azimuth is measured from
!                a fixed direction of the detector.
!   path: output. length from the surface of the earth to the detector
!                 in g/cm**2 for given tet, fai.
       subroutine muroc0(dsn, a, z, zba, z2ba, rho, bh, bv, beta)
         parameter (nzen=20, nazim=100)
         parameter (pi=3.141592, Torad=pi/180.)
         dimension dh(nzen, nazim), za(nzen), aa(nazim)
!
             character*(*) dsn
             data jflat/0/
             character*79 txt
!
             open(21, file=dsn, status='shr', action='read')
             write(*,*) ' rock profile is by ', dsn
             read(21,'(a)') txt
             write(*,'(a)') txt
!
             read(21,'(a)') txt
             write(*,'(a)') txt
             read(21,'(a)') txt
             write(*,'(a)') txt
             read(txt,*)  a, z, zba, z2ba, rho, bh, bv, beta
             rhoi=rho
!
             read(21,'(a)') txt
             l=index(txt, 'flat')
             if(l .gt. 0) then
                 read(txt(l+5:20), *) vdep
                 write(*,*) ' flat earth surface assumed ',
     *           ' vertical depth=',vdep,' hg/cm**2'
                 jflat=1
                 vdep=vdep*100.
             else
                 jflat=2
                 call ksetrv(za, 1, nzen, -1.)
                 read(21,*) za
                 call kfrge(za, 1, -nzen, 0., nzena, icon)
                 nazima=0
                  do   i=1, nazim
                     read(21, *, end=101) aa(i), (dh(j, i),j=1, nzena)
                      do   j=1, nzena
!                           convert dh into g/cm**2
                        dh(j,i)=dh(j,i)*100.*rhoi
                      enddo
                     nazima=nazima+1
                  enddo
  101            continue
                 dz=za(nzena)/(nzena-1)
                 da=360./nazima
                 write(*,*) ' zenith step=',dz, ' azimuth step=',da,
     *           ' # of zenith =',nzena, ' # of azimuth=',nazima
                 dz=dz*Torad
                 da=da*Torad
            endif
            close(21)
            return
!
       entry murocq(zmin, zmax)
            if(jflat .eq. 1) then
                zmin=vdep
                zmax=vdep/.5
            elseif(jflat .eq. 2) then
                 zmin=dh(1,1)
                 zmax=zmin
                  do   i=1, nazima
                     call kfmin(dh, 1, nzena, lmn)
                     call kfmax(dh, 1, nzena, lmx)
                     zmin=min(dh(lmn, i), zmin)
                     zmax=max(dh(lmx, i), zmax)
                  enddo
            endif
            return
!
       entry murock(tet, fai, path)
!
            if(jflat .eq. 2) then
!                     dh is now in g/cm**2
                call
     *          k4ptdi(dh, nzena, nazima, nzen,
     *          0., 0.,  dz, da,  tet, fai, path)
            elseif(jflat .eq. 1) then
                path=vdep/cos(tet)
            else
                write(*,*) ' muroc0 not yet called '
                stop
            endif
       end
!     ****************************************************************
!     *                                                              *
!     * muradl:  compute radiation length of given rock              *
!     *                                                              *
!     *********************** tested 90.03.01 ************************
!
!   /usage/
!            call kradl1(z, zba,  z2ba, x0g)
!
!    z:  <z> charge of the matter
!  zba:  <z/a>
! z2ba: <z**2/a>
!  x0g:  //                           g/cm**2
!
!
!     *** note ***
!
!         correction to born approximation is not included in this
!         r.l so it must be included in the cross-section.
!
!
!
      subroutine muradl(z, zba, z2ba, x0g)
!
!         cnst=
!         4/137* r0**2 * n  where r0 is the classical electron radius
!                           n the avogadro number
!                           r0=2.8176e-13 cm
!                           n=6.0247
!
      data  cnst/1.396e-3/
!
      z3=z**(-0.3333333)
      alogz3=log( 183.* z3 )
      gzai=log(1440. * z3**2 ) /  alogz3
!        inverse of r.l in g/cm**2
!     t0inv=cnst / a  *  z*( z + gzai ) * alogz3
      t0inv=cnst * (z2ba+zba*gzai) * alogz3
!
      x0g=1. / t0inv
      end
!      ***********************************************************
!      *
!      * kmover: move real*4 array to another array
!      *
!      ***********************************************************
!/usage/ call kmover(a, intv, n, b, intb)
!         a--->  b
!
       subroutine kmover(a, intv, n, b, intb)
!
          dimension a(intv, n), b(intb, n)
!
           do   i=1, n
              b(1, i)=a(1, i)
           enddo
       end
      subroutine kcossn(cs,sn)
!     to generate cos(phy),sin(phy), where phy is uniform in (0,2pi).
!     .........required subprogram...rndc ..............................
!      *** until loop*** 
      do while (.true.)
          call rndc(u)
          call rndc(v)
          v=v+v-1.
          a=u*u
          b=v*v
          c=a+b
      if         (c .le. 1.)
     *                   goto 10
      enddo
   10 continue
      cs=(a-b)/c
      sn=2.*u*v/c
      return
      end
!     ****************************************************************
!     *                                                              *
!     *   kgauss: generates 2 independent gaussian random variables  *
!     *                                                              *
!     ****************************************************************
!
!  /usage/
!                call kgauss(av, s, x1, x2)
!
!     av:  average of gaussian distribution
!      s:  standard deviation //
!  x1,x2:  obtained 2 independent gaussian random variables
!
!
!   method:
!             generate following two:
!             sqrt(-log(u)) * cos(p) and  sqrt(-log(u)) * sin(p)
!             then, they are independent gaussian radom variables
!             with mean 0 and variance 1.  here, u is uniform random
!             number in (0,1), and p in (0,2pi)
!        *** note ***  establish condition for rndc usage, if want to
!                      select one of rnd or rnd2.
!
!
      subroutine kgauss(av,s,x1,x2)
!
!
!         *** until loop*** 
         do while (.true.)
             call rndc(u1)
             call rndc(u2)
!              generate cos(p) and sin(p)
             u1=u1+u1-1.
             u1s=u1*u1
             u2s=u2*u2
             tmp=u1s+u2s
         if         (tmp .lt. 1.)
     *                      goto 20
         enddo
   20    continue
         call rndc(u3)
         al=alog(u3)
         al=sqrt(-al-al)
         cs=(u1s-u2s)/tmp
         sn=u1*u2/tmp
         sn=sn+sn
         x1=al*cs*s+av
         x2=al*sn*s+av
      end
!     ****************************************************************
!     *                                                              *
!     * k4ptdi:  4-point two dimensional interpolation               *
!     *                                                              *
!     ****************************************************************
!
!   /usage/
!
!       call
!              k4ptdi(f, im, jm, iadj,  x0, y0, hx, hy, x, y, ans)
!
!     f is a 2-dimensional table of some function with 2 arguments.
!     f containes the function values at (x0,y0), (x0+hx, y0+hy),...
!     (x0+(im-1)*hx, y0+(jm-1)*hy).
!     iadj is the adjustable dimension.
!     ans gets the value of the funtion at (x,y)
!
!
!
      subroutine k4ptdi(f, im, jm, iadj, x0, y0, hx, hy, x, y, ans)
      dimension f(iadj,jm)
!
      a=(x-x0)/hx
      b=(y-y0)/hy
      i=a
      j=b
      i=min0(max0(i,0)+1,im-1)
      j=min0(max0(j,0)+1,jm-1)
      p=a+1.-i
      q=b+1.-j
      p1=1.-p
      q1=1.-q
      ans=( f(i,j)*p1 + f(i+1,j)*p ) * q1 +
     *                  ( f(i,j+1)*p1 + f(i+1,j+1)*p ) * q
      return
      end
!             get theta and fai of direction cos. in rad
       subroutine mudtoa(vx, vy, vz, teta, fai)
           parameter (pi=3.141592, Todeg=180./pi)
           real*8 vx, vy, vz
           if(vz .gt. 1.d0) then
              teta=0.
           else
              teta=acos(vz)
           endif
           if(abs(teta) .lt. 1.e-4) then
               fai=0.
           else
               fai=atan2(vy, vx)
           endif
       end
