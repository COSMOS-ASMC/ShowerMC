       include  'flux0.f'
!      ***************************************************************
!      *
!      * mudet:  muon bundle observation : process the result of
!      *         mufug assuming a paticular detector structure.
!      *
!      *  /usage/   prepare data which is output from mufug.
!      *         the only one module to be imported is mumin which
!      *         may be in #load.load if you have been using mufug.
!      *         in that case, fort command may be used to run the job.
!      *       ( mylib lib1(#load.load)
!      *         alloc f(sysinc) da('c2g5100.cosmos.gem') shr
!      *         is  needed at the session start time).
!      *   note:  input and output for input and output dataset
!      *          names cannot be changed at run time.  these
!      *          must be chagned below in the data statement.
!      *
!      *   if mufug has not been used, you may add *include statement
!      *   before the first line of the program to include mumin.
!      *   the parameters can be put through namelist at run time or
!      *   you may change block data in the last part of this file.
!      *   (if want to use defaults, put /* for read action).
!      * --------------------------------------------------------
!      * for the meaning of the parameters, read block data.
!      * --------------------------------------------------------
!      *
!      *         this program may be run at tss mode because the
!      *         run time is not so long.
!      *
!      *  the output from this program is special for gd command
!      * (heading, and caption for 5 items exist but no stop code).
!      *  5 items written in each record are:
!      *
!      *       nmu, e1ry, heavy index, cos(1ry), z1
!      * where
!      *    nmu: # of muons in an event
!      *   e1ry: 1ry energy in tev
!      * heavy index: 1-7 for p,alfa,l,m,h,vh,fe
!      *  cos(1ry):   cos of zenith angle of 1ry
!      *    z1:       first collision point of 1ry (g/cm**2)
!      *
!      *
!      *    you may want to have different output (say, enegy
!      *    of each muon,  distance of pair of muons etc).
!      *    in that case, you may make another program modifing
!      *    this one.
!      *
!      *  to export this one to other machine,
!      *     $mudet
!      *     $flux
!      *     mudet
!      *     mumin
!      *  are needed.
!      *
!      ***************************************************************
!             input file
        character*24 input/'c2s5001.#gd2.data'/
!             output file
        character*24 output/'c2s5001.#gd.data'/
!
!                 input data file
        open(7, file=input, action='read',
     *           status='shr')
!                 output file
        open(8,  file=output, action='write')
!          init
       call mudet0
!          execute detection
       call mudetx
!          end of all work
       call mudete
      end
      subroutine mudet0
       include  'Zflux.f'
       include  'Zmudet.f'
        character*90 txt
        character*16 cap(5)
!       data cap/'n muon',' e 1ry(tev)', '1ry indx', 'cos(zenith)',
        data cap/'n muon',' e 1ry(tev)', '1ry indx', 'zen angl(deg)',
     *           'z1(g/cm**2)'/
        call muprm(6,'w', icon)
!            read parameter
        call muprm(5, 'r', icon)
!              get minimum energy to reach the depth
        call mumin(depth, ec)
!              size of area to let the 1ry fall (cm)
        fallx=detx + 40.e2/ ec
        fally=dety + 40.e2/ ec
        if(jcent .eq. 0 .and. n1ry .ne. 1) then
            write(*,*) ' n1ry is made to be 1 because jcent=0'
            n1ry=1
        endif
!         init resample procedure
!              dco: coeff. for original dif. flux
!              dcr: //         reampled //
        call  flux0(s1,s2,s1r,s2r,s12r, scut, dco, dcr)
!           write parameter
        call muprm(6,'w',icon)
        write(*,*) ' minimum energy to reach the det.=',ec, ' tev'
        write(txt,'('' id='',a,'' site='',a, ''area of 1ry='',f6.2,
     *   '' m x'',  f6.2, '' m.  n1ry='',i5, '' jcent='',i1)')
     *  model, site, fallx/100., fally/100., n1ry, jcent
        write(8) txt
        write(8) cap
        mucmax=0
      end
!
      subroutine mudetx
       include  'Zflux.f'
       include  'Zmudet.f'
!
         do   i=1, 9999999
           read(7, end=200)
     *     nshwno, e1ry, ihg1ry, wx, wy, wz,zfirst
!              resampling and/or e0cut
           call mursmp(jcon)
           if(jcon .ne. 0) then
!               rejected, discard event
              read(7)
           else
              read(7)
     *        muc, (eats(j), ed(j), xd(j), yd(j), txd(j), tyd(j),
     *        tzd(j),   j=1, muc)
!             do 5 kkk=1, muc
!                  write(14) eats(kkk), ed(kkk)
!   5         continue
!                see if accross the detector
              call mudeti
           endif
         enddo
  200   continue
      end
      subroutine mudeti
!          observation by square detector
       include  'Zflux.f'
       include  'Zmudet.f'
        mucmax=max(muc, mucmax)
!          get position of ptcl at z=detz and z=0.
         do   i=1, muc
!             position at z=detz
            xt(i)=xd(i)
            yt(i)=yd(i)
!              position at z=0.
            t=- detz/tzd(i)
            xb(i)=xd(i) + t* txd(i)
            yb(i)=yd(i) + t* tyd(i)
         enddo
!           this event is equivalent to n1ry primaries
         do   n=1, n1ry
!             bndle counter
           nbndl=0
!            falling place of 1ry
           call mu1ryl(n, xc, yc)
!              every ptcl position should be displaced by (xc,yc)
!              detector should be move to (0,0,0) which was assumed
!              at -(detx,dety,detz)/2.  then xc+detx/2, etc must be
!              added
           dx=xc+detx/2
           dy=yc+dety/2
            do   i=1, muc
               xtt=xt(i)+dx
               ytt=yt(i)+dy
               xbt=xb(i)+dx
               ybt=yb(i)+dy
!                 see if the segment accross the detector
               call k3dclp(xtt,ytt,detz, xbt, ybt,0., detx, dety,
     *         detz, x1,y1,z1, x2,y2,z2, icon)
               if(icon .eq. 0) then
!                  accross; compute length
                  path=sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
                  if(path .gt. pathmn) then
                      nbndl=nbndl+1
                  endif
               endif
            enddo
           if(nbndl .gt. 0) then
!               output
!             write(8) float(nbndl), e1ry, float(ihg1ry), wz, zfirst
              tet=acos(wz)*180/3.1415
              write(8) float(nbndl), e1ry, float(ihg1ry),tet, zfirst
           endif
         enddo
       end
       subroutine mudete
       include  'Zflux.f'
       include  'Zmudet.f'
          write(*,*) ' end of all events '
          write(*,*) ' max mu #=',mucmax
       end
!
       subroutine mu1ryl(n, xc,yc)
       include  'Zflux.f'
       include  'Zmudet.f'
!           samples 1ry location
           if(jcent .eq. 0) then
               xc=0.
               yc=0.
           else
               call rndc(ux)
               call rndc(uy)
               if(mod(n,2) .eq. 0) then
                  xc=(2*ux-1.)*fallx
                  yc=(2*uy-1.)*fally/wz
               else
                  xc=(2*ux-1.)*fallx/wz
                  yc=(2*uy-1.)*fally
               endif
           endif
       end
!       resampling routine
       subroutine mursmp(jcon)
!               resample 1ry energy by rejection
!           or cut 1ry > e0cut*mass
       include  'Zflux.f'
       include  'Zmudet.f'
        if(resmpl) then
             call ar1rye(dco(ihg1ry), beta1(ihg1ry), beta2(ihg1ry),
     *       ebenda(ihg1ry), gbend,
     *       dcr(ihg1ry),  abeta1(ihg1ry), abeta2(ihg1ry),
     *       aebend(ihg1ry), agbend,  e1ry, jcon)
        else
            jcon=0
        endif
        if(e1ry .gt. e0cut) then
            jcon=1
        endif
       end
!     ******************************************************************
!     *                                                                *
!     *  ar1rye: accept or reject 1ry energy                           *
!     *                                                                *
!     *********************** tested 83.03.01 **************************
!
!        /usage/  call ar1rye( dco, b1, b2, eb1, gb1,
!                              dcr, g1, g2, eb2, gb2, e, icon)
!
!      given an energy e which has been sampled according to
!      a 1ry spectrum speicfied by (b1,b2,eb1,gb1).
!      this subroutine determines that e be accepted or rejected
!      so that accepted ones become a 1ry specified by
!      (g1,g2,eb2,gb2).
!
!  -- input --
!     dco:  differential coeff. to be used for original 1ry
!           (dco and dcr can be obtained by calling flux0)
!      b1:  slope of 1ry at energy < eb1
!      b2:  //                     >
!     eb1:  bending point of 1ry
!     gb1:  logical.
!           if gb1=t then
!
!               f(e)de= e**(-1-b1) * (1 + e/eb)**(b1-b2) de
!
!               is assumed.    slope at e=eb is just -(b1+b2)/2.
!
!           else
!                      !    e**(-1-b1)de                   e<eb1
!               f(e)de=!
!                      !    e**(-1-b2) eb**(b2-b1) de      e>eb1
!
!               is assumed.
!
!
!    dcr: same as dco for resampled one.
!   g1,g2:  same as b1,b2.  for object 1ry
!     eb2:     same as eb1.    //
!     gb2:     same as gb1.    //
!
!       e:  energy
!
!  *** note ***
!      if b1=b2, eb1 and gb1 are not used. f(e)de=e**(-1-b1) de
!      is assumed.
!
!
! -- output --
!    icon:   if 0, e is to be accepted.  if 1, e is to be rejected.
!
!  *** note ***
!        rnd1 or rnd2 must be fixed beforehand
!        in the case of gb2=t,  it is desirable to be g1>=b1 else
!        result is not so much accurate.
!
!
      subroutine ar1rye(dco, b1, b2, eb1, gb1,
     *                  dcr, g1, g2, eb2, gb2, e, icon)
!
!
      logical gb1, gb2
!
!
           a=fd1ry(b1,b2,eb1,gb1,e)*dco
           b=fd1ry(g1,g2,eb2,gb2,e)*dcr
!               flux ratio at e  (2nd/1st)
           r=b/a
           call rndc(u)
           if(u .lt. r) then
               icon=0
           else
               icon=1
           endif
      end
!     ******************************************************************
!     *                                                                *
!     * fd1ry: gives differntial 1ry flux                              *
!     *                                                                *
!     ***********************  tested 83.02.28 *************************
!
!    /usage/   flux=fd1ry(b1,b2,eb, gbend, e)
!
!  -- input --
!     b1:  1ry index at e<eb
!     b2:  //            >
!     eb:  bending point of 1ry
!  gbend:  logical.   to show gradually bending 1ry or not.
!
!          if gbend=f  then
!
!                      !  e**-(b1+1) de     at e<eb
!               f(e)de=!
!                      ! eb**(b2-b1) * e**-(b2+1) de     at e>eb
!
!          else
!
!               f(e)de= e**-(b1+1) * (1+e/eb)**(b1-b2) de
!
!          is assumed
!
!    *** note ***
!          if b1 = b2  then eb and gbend are not used and single slope
!          1ry is assumed:  f(e)de= e**-(b1+1)  de
!
!
!  -- process --
!          gives f(e)
!
!  -- output --
!       this is a function program.  unit of flux is arbitray.
!
!                              =   =   =   =
!
                    function fd1ry( b1, b2, eb, gbend, e )
!
      logical gbend
!
      if( b1 .eq. b2 ) then
          fd1ry=e**(-b1-1.)
      elseif(gbend) then
          fd1ry=e**(-b1-1.) * (1. + e/eb)**(b1 - b2)
      elseif( e .lt. eb) then
          fd1ry=e**(-b1-1.)
      else
          fd1ry=eb**(b2-b1) * e**(-b2-1.)
      endif
      return
      end
      subroutine muprm(io, rw, icon)
       include  'Zflux.f'
       include  'Zmudet.f'
         character*(*) rw
         namelist /muparm/
     * abeta1, abeta2, aebend, agbend, beta1, beta2, depth, detx,
     * dety, detz, ebenda, enorm, e0cut, e1, fallx, fally, gbend,
     * n1ry, pathmn, resmpl, site, model

!
       if(rw .eq. 'w') then
          write(io, muparm)
          icon=0
       elseif(rw .eq. 'r') then
          read(io, muparm, end=800)
          icon=0
       else
          stop ' bad rw for muprm '
       endif
       e2=e0cut
       return
  800  continue
       icon=1
       end
!         test k3dclp
!     data a/1./, b/1./, c/1./
!     x0=.5
!     y0=-.5
!     z0=1.
!     x1=.0
!     y1=1.5
!     z1=0.0
!          call k3dclp(x0, y0, z0, x1, y1, z1, a, b, c,
!    *            xo0, yo0, zo0, xo1, yo1, zo1, icon)
!      write(*,*) ' icon=',icon, ' xo0=',xo0, ' yo0=',yo0, ' zo0=',zo0
!      write(*,*) '        ', ' xo1=',xo1, ' yo1=',yo1, ' zo1=',zo1
!      end
!        *********************************************************
!        *
!        * k3dclp: 3d clipping of a line segment by a box
!        *
!        *********************************************************
!
! /usage/  call k3dclp(xi0, yi0, zi0, xi1, yi1, zi1, a, b, c,
!         *            xo0, yo0, zo0, xo1, yo1, zo1, icon)
!
!       xi0, yi0, zi0: input.  1st point of the line segment
!       xi1, yi1, zi1: input.  2nd point of the line segment
!         a, b, c: input.  the lenght of three edges of the box.
!                   the three edges lie on (0,0,0)-(a,0,0)
!                                          (0,0,0)-(0,b,0)
!                                          (0,0,0)-(0,0,c)
!       xo0,yo0,zo0: output.  1st point of the segment
!       xo1,yo1,zo1: output.  2nd point of the segemnt
!                             these may be orignal points
!                             if points are contained in the
!                             box, or the point(s) on the surface
!                             of the box where the segment accrosses
!                             the box or not given.
!       icon: ouput.          =-1 --> no crossing point at all.
!                             = 0 --> crossing points exist or contained
!
       subroutine k3dclp(xi0, yi0, zi0, xi1, yi1, zi1, a, b, c,
     *                   xo0, yo0, zo0, xo1, yo1, zo1, icon)
!
!               ix : to set bit at x-th position  x=1,2, from right
!
       parameter (i6=5, i5=4, i4=3, i3=2, i2=1, i1=0)
       integer bp6, bp5, bp4, bp3, bp2, bp1
!              bit pattern whose x-th bit is on, where x is bpx)
       parameter (bp6=2**i6, bp5=2**i5, bp4=2**i4,
     *            bp3=2**i3, bp2=2**i2, bp1=2**i1)
       logical ok
!
       jc0=0
       jc1=0
!
       x0=xi0
       y0=yi0
       z0=zi0
       x1=xi1
       y1=yi1
       z1=zi1
!
       if(z0 .gt. c ) then
           jc0=ibset(jc0, i6)
       endif
       if(z0 .lt. 0.) then
           jc0=ibset(jc0, i5)
       endif
       if(y0 .gt. b) then
           jc0=ibset(jc0, i4)
       endif
       if(y0 .lt. 0.) then
           jc0=ibset(jc0, i3)
       endif
       if(x0 .gt. a) then
           jc0=ibset(jc0, i2)
       endif
       if(x0 .lt. 0.) then
           jc0=ibset(jc0, i1)
       endif
!
       if(z1 .gt. c ) then
           jc1=ibset(jc1, i6)
       endif
       if(z1 .lt. 0.) then
           jc1=ibset(jc1, i5)
       endif
       if(y1 .gt. b) then
           jc1=ibset(jc1, i4)
       endif
       if(y1 .lt. 0.) then
           jc1=ibset(jc1, i3)
       endif
       if(x1 .gt. a) then
           jc1=ibset(jc1, i2)
       endif
       if(x1 .lt. 0.) then
           jc1=ibset(jc1, i1)
       endif
!
!       *** until loop*** 
       do while (.true.)
           ok=jc0 .eq. 0 .and. jc1 .eq. 0
           if(ok) then
               icon=0
               xo0=x0
               yo0=y0
               zo0=z0
               xo1=x1
               yo1=y1
               zo1=z1
           else
               ok=iand(jc0, jc1) .ne. 0
               if(ok) then
                   icon=-1
               else
                   if(jc0 .eq. 0) then
                      tmpx=x0
                      tmpy=y0
                      tmpz=z0
                      x0=x1
                      y0=y1
                      z0=z1
                      x1=tmpx
                      y1=tmpy
                      z1=tmpz
                      itmp=jc0
                      jc0=jc1
                      jc1=itmp
                   endif
                   if(iand(jc0, bp6) .ne. 0) then
                       t=(c-z0)/(z1-z0)
                       z0=c
                       x0=x0 + (x1-x0) * t
                       y0=y0 + (y1-y0) * t
                   elseif(iand(jc0,bp5) .ne. 0) then
                       t=  -z0/(z1-z0)
                       z0=0.
                       x0=x0 + (x1-x0) * t
                       y0=y0 + (y1-y0) * t
                   elseif(iand(jc0, bp4) .ne. 0) then
                       t=(b-y0)/(y1-y0)
                       y0=b
                       x0=x0 + (x1-x0) * t
                       z0=z0 + (z1-z0) * t
                   elseif(iand(jc0, bp3) .ne. 0) then
                       t=  -y0/(y1-y0)
                       y0=0.
                       x0=x0 + (x1-x0) * t
                       z0=z0 + (z1-z0) * t
                   elseif(iand(jc0, bp2) .ne. 0) then
                       t=(a-x0)/(x1-x0)
                       x0=a
                       y0=y0 + (y1-y0) * t
                       z0=z0 + (z1-z0) * t
                   else
                       t=  -x0/(x1-x0)
                       x0=0.
                       y0=y0 + (y1-y0) * t
                       z0=z0 + (z1-z0) * t
                   endif
!
                   jc0=0
                   if(z0 .gt. c ) then
                       jc0=ibset(jc0, i6)
                   endif
                   if(z0 .lt. 0.) then
                       jc0=ibset(jc0, i5)
                   endif
                   if(y0 .gt. b) then
                       jc0=ibset(jc0, i4)
                   endif
                   if(y0 .lt. 0.) then
                       jc0=ibset(jc0, i3)
                   endif
                   if(x0 .gt. a) then
                       jc0=ibset(jc0, i2)
                   endif
                   if(x0 .lt. 0.) then
                       jc0=ibset(jc0, i1)
                   endif
               endif
           endif
       if         (ok)
     *                    goto 100
       enddo
  100  continue
      end
       block data
       include  'Zflux.f'
       include  'Zmudet.f'
!          following must be fixed here or at run time by namelist
!       ---------------------------------------------------------------
!                     original spectrum
!                 p      alfa     l      m(cno)   h       vh     fe
           data
     1     beta1/ 1.7,    1.7,   1.7,   1.6,     1.6,     1.6,   1.4/,
     2     beta2/ 2.,     2.,    2.,    2.,      2.,      2.,     2./,
     3     ebenda/3000., 3000., 3000., 3000., 3000.,  3000.,   3000./
           data  gbend/.true./
!           aflux: integral flux value (in some unit)
!           eflux: total kinetic energy where aflux is given (tev)
      data aflux/ 42.,    20.,   .6,    14.,     10.,     4.,    10./,
     1     eflux/ .56,    .56,   .56,   .56,     .56,     .56,   .56/
!
!
!                     resamling spectrum
           data
     1    abeta1/ 1.7,    1.7,   1.7,   1.6,     1.6,     1.6,   1.5/,
     2    abeta2/ 2.,     2.,    2.,    2.,      2.,      2.,     2./,
     3    aebend/ 100.,  200.,  400.,  700., 1200.,  1700.,   2600./,
     4    enorm/.56, .56, .56, .56, .56, .56, .56/
!    1    abeta1/ 1.7,    1.7,   1.7,   1.7,     1.7,     1.7,   1.7/,
!    3    aebend/ 300.,  600., 1200., 2100., 3600.,  5100.,   7800./,
!    3    aebend/ 7*3000./,
!             enorm is the minimum energy where the differential
!             flux value is made to be the same as the original one
!
          data agbend/.true./
!              resample or not    # of 1ry made from 1 event
          data resmpl/.true. /, n1ry/ 500/
!             for identifing the site later (comment) (6 cha)
          data site/'frj   '/, model/'p-poor(1)'/
!             detector edge size (frj)
          data detx/1230./, dety/600./,  detz/ 600./
!             detector edge size (kgf)
!         data detx/ 600./, dety/650./,  detz/ 600./
!               minimum path length for detection (cm)
          data pathmn/100./
!           approximate  minimum  detector depth (  in g/cm**2 )
          data depth/4500.e2/
!           if a higer energy region is covered by different
!           set of simulation which is done by einp=e0cut,
!           give that value (heavy is cut at a*e0cut)
          data e0cut/ 1000000./
!           e1 = emin of 1ry ( e2=e0cut is set in the program)
          data e1/  3000./
!            if want to fall the 1ry in the center of the detector
!            make jcent=0. in this case, n1ry is made to 1.
!            if detx and dety is large, and pathmn= small, you can
!            see the original muon distribution at the site.
          data jcent/1/
!      ---------------------------------------------------------------
       end
