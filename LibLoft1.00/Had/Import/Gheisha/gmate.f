c      subroutine gsmate(imat,namate,a,z,dens,radl,absl)
      subroutine gsmate(imat,a,z,dens,radl,absl)
      implicit none
c       integer imat: input.  matter
c       character namate: input  matter name
c       real      a:input.  mass number
c       real      z:input.  charge
c       real   dens:input.  density (g/cm3) 
c       real   radl:input.  radiaiton length  (cm)
c       real   absl:input. 
      integer *4 maxmat,maxmix
      parameter (maxmat=50)
      parameter (maxmix=5)
      character *20 name(maxmat),nam20
      real *4 amat(maxmat)
      real *4 zmat(maxmat)
      real *4 densmat(maxmat)
      real *4 abslmat(maxmat)
      real*4 aheffmat(maxmat)
      integer *4 nmixmat(maxmat)
      integer *4 imixmat(maxmix,maxmat)
      real *4 wmixmat(maxmix,maxmat)
      real *4 amix(maxmix)
c
      integer *4 imat,flag,nmix,imix(1)
      real *4 a,z,dens,absl,radl,wmix(1)
c     character *(*) namate
      character *10 c10
      integer nlmat
      real*4 wmat(1)
      integer *4 what(1)
      integer *4 i,j,start,finish,nlm
      real *4 ghsigm,abseff
      real *4 amol,zmol,aeff,zeff,aheff
      integer *4 fail
c
      save name
      save amat,zmat,densmat,abslmat,aheffmat,nmixmat,imixmat
      save wmixmat
      data amat/maxmat*0.0/
      data j/0/
      data nlm/0/
c
      fail=1301
      if(imat.le.0) go to 1313
      if(imat.gt.maxmat) go to 1313
c     name(imat)=namate
      amat(imat)=a
      zmat(imat)=z
      densmat(imat)=dens
      abslmat(imat)=absl
      nmixmat(imat)=1
      imixmat(1,imat)=imat !contains only itself
      aheffmat(imat)=a
c
      return
c
      entry gfmate(imat,nam20,a,z,dens,radl,absl,nmix,imix,wmix)
c
      fail=1302
      if(imat.le.0) go to 1313
      if(imat.gt.maxmat) go to 1313
      nam20=name(imat)
      a=amat(imat)
      z=zmat(imat)
      dens=densmat(imat)
      radl=1.   !no need for rad. lengths here
      absl=abslmat(imat)
      fail=1308
      if(nmix.lt.nmixmat(imat)) go to 1313
      nmix=nmixmat(imat)
      call ucopy(imixmat(1,imat),imix,nmix)
      call ucopy(wmixmat(1,imat),wmix,nmix)
c
      return
c
c     entry gsmixt(imat,namate,nlmat,wmat,what,dens)
      entry gsmixt(imat,nlmat,wmat,what,dens)
c.
c.    ******************************************************************
c.    *                                                                *
c.    *       defines mixture or compound imat as composed by          *
c.    *       the basic nlmat materials defined by arrays a,z and wmat *
c.    *                                                                *
c.    *       if nlmat.gt.0 then wmat contains the proportion by       *
c.    *       weigths of each basic material in the mixture.           *
c.    *                                                                *
c.    *       if nlmat.lt.0 then wmat contains the number of atoms     *
c.    *       of a given kind into the molecule of the compound        *
c.    *       in this case, wmat in output is changed to relative      *
c.    *       weigths.  not used                                       *
c.    *                                                                *
c.    *       nb : the radiation length is computed according          *
c.    *            the egs manual slac-210 uc-32 june-78               *
c.    *                           formula  2-6-8 (37)                  *
c.    *                                                                *
c.    *    ==>called by : <user>, ugeom                                *
c.    *       authors    r.brun, m.maire  *********                    *
c.    *                                                                *
c.    ******************************************************************
c.
c.
c.    ------------------------------------------------------------------
c.
      fail=1303
      if (imat.le.0) go to 1313
      if(imat.gt.maxmat) go to 1313
c
c             compute proportion by weigths in the compound
c
      nlm=iabs(nlmat)
      fail=1306
      if(nlm.le.1.or.nlm.gt.maxmix) go to 1313
      call ucopy(what,imixmat(1,imat),nlm)
      if(nlmat.lt.0) then
         amol   = 0.
         zmol   = 0.
         do 10 i= 1,nlm
          j=what(i)
          fail=1304
          if(j.lt.0.or.j.gt.maxmat) go to 1313
          fail=1305
          if(amat(j).eq.0.0) go to 1313
          amol   = amol + wmat(i)*amat(j)
          zmol   = zmol + wmat(i)*zmat(j)
   10    continue
         do 20 i= 1,nlm
          j=what(i)
          wmat(i)= wmat(i)*amat(j) / amol
   20    continue
      endif
c
c             compute effective mixture parameters
c
      aeff   = 0.
      zeff   = 0.
      do 40 i = 1,nlm
         j=what(i)
         fail=1304
         if(j.lt.0.or.j.gt.maxmat) go to 1313
         fail=1304
         if(abslmat(j).eq.0.0) go to 1313
         amix(i)=amat(j)
         aeff   = aeff + wmat(i)*amat(j)
         zeff   = zeff + wmat(i)*zmat(j)
   40 continue
      call ghmix(amix,wmat,nlm,aheff)
      abseff=10000.*aheff/(6.022*dens*ghsigm(5.,5,aheff))
c
c     name(imat)=namate
      amat(imat)=aeff
      zmat(imat)=zeff
      densmat(imat)=dens
      abslmat(imat)=abseff
      nmixmat(imat)=nlm
      aheffmat(imat)=aheff
      call ucopy(wmat,wmixmat(1,imat),nlm)
c
      return
c
      entry gmatep(flag)
       fail=1307
       if(flag.lt.0.or.flag.gt.maxmat) go to 1313
       start=1
       finish=maxmat
       if(flag.ge.1) then
        start=flag
        finish=flag
       endif !
       write(6,1002)
1002   format
     x(1x,'mat name                        a          z        rl' )
       do 1000 i=start,finish
        if(amat(i).eq.0.0) go to 1000
        write(6,1001) i,name(i),amat(i),zmat(i),c10
1001    format(1x,i3,1x,a20,2f10.3,a10)
1000   continue
      return
c
1313  continue
      write(6,*) 'fail=',fail,' in gsmate'
      if(fail.eq.1301) then
       write(6,*) 'imat=',imat,' is out of range'
       write(6,*) 'this material remains undefined'
      else if(fail.eq.1302) then
       write(6,*) 'imat=',imat,' is out of range'
       write(6,*) 'cannot find material properties'
      else if(fail.eq.1303) then
       write(6,*) 'imat=',imat,'  is out of range'
       write(6,*) 'this material mix remains undefined'
      else if(fail.eq.1304) then
       write(6,*) 'imat=',j,' is out of range'
       write(6,*) 'this material cannot be used in a mixture'
      else if(fail.eq.1305) then
       write(6,*) 'imat=',j,' has not been defined'
       write(6,*) 'this material cannot be used in a mixture'
      else if(fail.eq.1306) then
       write(6,*) 'nlm=',nlm,' too many or too few mixture components'
       write(6,*) 'this material mix remains undefined'
      else if(fail.eq.1307) then
       write(6,*) 'flag=',flag,' is out of range'
       write(6,*) 'can not print'
      else if(fail.eq.1308) then
       write(6,*) 'nmix=',nmix,' not enough room for requested mixture'
       write(6,*) 'no mixture is returned'
      endif !fail.eq.1301
      return
      end
c
      subroutine ghmix(a, w, n, aeff)
c
c     ******************************************************************
c     *                                                                *
c     *   works out an effective atomic weight aeff for a material     *
c     * with n elements of atomic weight a(i) in proportion w(i) by    *
c     * weight.  the criterion is that the hadronic interaction        *
c     * length of a 5 gev/c pion is correct.  errors on the calculated *
c     * hadronic interaction length for other momenta and other        *
c     * particles in geant version 3.04 are less than 1% in most cases.*
c     * for details see memorandum opal/0037n/ja/md, ref. hadron       *
c     * milestone 84/003, calculation of hadronic interaction lengths  *
c     * for mixtures.                                                  *
c.    *                                                                *
c.    *    ==>called by : gsmixt                                       *
c.    *       author    j.allison  *********                           *
c     *                                                                *
c     ******************************************************************
c
c
      integer *4 n,fail
      real a(n), w(n),aeff
c.
c.    ------------------------------------------------------------------
c.
c         work out pint which is proportional to the interaction
c         probability.  also work out the mean atomic weight, i.e. that
c         weighted by proportion by numbers, as starting value for
c         iterative method of finding aeff.
c
      pint   = 0.
      ainv   = 0.
      wtot   = 0.
      do 10 i = 1, n
         pint   = pint + w(i) * ghsigm(5., 5, a(i)) / a(i)
         ainv   = ainv + w(i) / a(i)
         wtot   = wtot + w(i)
  10  continue
c
      fail=1301
      if ( abs ( wtot - 1. ) .gt. 0.01 ) go to 1313
c
c         work out aeff which gives pint for 5 gev/c pion.
c         (this is a short newton's method loop.)
c
      aeff   = 1. / ainv
      pnew   = ghsigm ( 5., 5, aeff ) / aeff
      da     = 1.
      iter   = 0
   20 continue
      iter   = iter + 1
      aeff   = aeff + da
      pold   = pnew
      daold  = da
      pnew   = ghsigm ( 5., 5, aeff ) / aeff
      dp     = pnew - pold
      da     = (pint - pnew ) * daold / dp
      fail=1302
      if ( ( iter .gt. 1 ) .and. ( abs ( da ) .gt. abs ( daold ) ) )
     +                            go to 1313
      if ( abs ( da ) .gt. 0.01 ) go to 20
c
      return
c
c            error conditions.
c
1313  continue
      write(6,*) 'fail=',fail,' in ghmix'
      if(fail.eq.1301) then
       write(6,*) 'mixtures weights do not add up to 1'
      else if(fail.eq.1302) then
       write(6,*) 'hadronic interaction mixture does not converge'
      endif
      end


      function ghsigm(p,ipart,a)
c.
c.    ******************************************************************
c.    *                                                                *
c.    *   returns absorption cross-section in millibarns               *
c.    *   for a particle with                                          *
c.    *                momentum p (gev/c)                              *
c.    *                gheish type ipart                               *
c.    *      on a nucleus of atomic number a                           *
c.    *                                                                *
c.    *       author    r.barlow  *********                            *
c.    *       modified to deal with k0 and extend to 41 particle       *
c.    *                types by john allison, 31/1/84. *               *
c.    *       modified to remove geant ref's by w.b.atwood             *
c.    *                15-apr-1993                                     *
c.    *                                                                *
c.    ******************************************************************
c.
c.             the k0 is internal type 6 which causes the
c              mean of the k+ and k- cross sections to be calculated.
c
c.
c.    ------------------------------------------------------------------
c.
      ghsigm = 1.e-20
      if (ipart .le. 0) return
      if (ipart .ge. 7) return
c
c      if (itype .ne. 6)then
      if (ipart .ne. 6)then
         ghsigm = ghsig(p, ipart, a)
      else
         ghsigm = 0.5 * (ghsig(p, 3, a) + ghsig(p, 4, a))
      endif
c
      end
c
      function ghsig(p,itype,a)
c.
c.    ******************************************************************
c.    *                                                                *
c.    *   returns absorption cross-section in millibarns               *
c.    *   for a particle with                                          *
c.    *                momentum p (gev/c)                              *
c.    *                type itype                                      *
c.    *      on a nucleus of atomic number a                           *
c.    *                                                                *
c.    *    ==>called by : ghsigm                                       *
c.    *       author    r.barlow  *********                            *
c.    *                                                                *
c.    ******************************************************************
c
c              the internal particle types are as follows...
c                1 proton
c                2 antiproton
c                3 k+
c                4 k-
c                5 pion
c
      dimension sigb(100),alpha(100),sig(162,5)
      dimension sigpr(162,1),sigpb(162,1),sigkp(162,1),sigkm(162,1)
      dimension sigpi(162,1)
      equivalence (sigpr(1,1),sig(1,1)),(sigpb(1,1),sig(1,2))
      equivalence (sigkp(1,1),sig(1,3)),(sigkm(1,1),sig(1,4))
      equivalence (sigpi(1,1),sig(1,5))
c
           data alpha/
     +  0.9826,0.9659,0.9500,0.9348,0.9203,0.9064,0.8932,0.8807,0.8687,
     +  0.8574,0.8466,0.8363,0.8265,0.8172,0.8084,0.8000,0.7921,0.7845,
     +  0.7773,0.7705,0.7640,0.7579,0.7520,0.7464,0.7411,0.7361,0.7313,
     +  0.7267,0.7224,0.7182,0.7142,0.7105,0.7069,0.7034,0.7001,0.6970,
     +  0.6939,0.6911,0.6883,0.6856,0.6831,0.6807,0.6783,0.6761,0.6739,
     +  0.6718,0.6698,0.6679,0.6660,0.6642,0.6625,0.6609,0.6592,0.6577,
     +  0.6562,0.6547,0.6533,0.6519,0.6506,0.6493,0.6481,0.6469,0.6457,
     +  0.6446,0.6435,0.6424,0.6413,0.6403,0.6393,0.6384,0.6374,0.6365,
     +  0.6356,0.6347,0.6339,0.6330,0.6322,0.6314,0.6307,0.6299,0.6292,
     +  0.6284,0.6277,0.6270,0.6263,0.6257,0.6250,0.6244,0.6237,0.6231,
     +  0.6225,0.6219,0.6213,0.6208,0.6202,0.6197,0.6191,0.6186,0.6181,
     +  0.6175/
           data sigb/
     +    1.02,  2.10,  3.21,  4.36,  5.54,  6.75,  7.99,  9.24, 10.51,
     +   11.79, 13.08, 14.37, 15.66, 16.95, 18.22, 19.49, 20.75, 21.99,
     +   23.22, 24.43, 25.62, 26.79, 27.93, 29.06, 30.16, 31.24, 32.30,
     +   33.33, 34.34, 35.32, 36.29, 37.22, 38.14, 39.03, 39.91, 40.75,
     +   41.58, 42.39, 43.18, 43.95, 44.70, 45.43, 46.14, 46.83, 47.51,
     +   48.17, 48.82, 49.45, 50.06, 50.66, 51.25, 51.83, 52.39, 52.93,
     +   53.47, 53.99, 54.51, 55.01, 55.50, 55.98, 56.45, 56.92, 57.37,
     +   57.81, 58.25, 58.68, 59.10, 59.51, 59.91, 60.31, 60.70, 61.08,
     +   61.46, 61.83, 62.19, 62.55, 62.90, 63.25, 63.59, 63.93, 64.26,
     +   64.58, 64.90, 65.22, 65.53, 65.84, 66.14, 66.44, 66.74, 67.03,
     +   67.31, 67.60, 67.88, 68.15, 68.42, 68.69, 68.96, 69.22, 69.48,
     +   69.74/
       data sigpr /
     +  79.01, 76.05, 73.81, 71.75, 69.67, 66.70, 64.89, 63.09, 60.58,
     +  58.99, 57.49, 56.09, 54.02, 52.75, 51.54, 50.39, 48.78, 47.79,
     +  46.86, 45.99, 44.73, 43.97, 43.25, 42.59, 41.59, 40.99, 40.41,
     +  39.85, 38.88, 38.26, 37.67, 37.06, 36.31, 35.82, 35.37, 35.00,
     +  34.40, 34.03, 33.71, 33.45, 33.00, 32.72, 32.50, 32.31, 32.01,
     +  31.84, 31.69, 31.56, 31.38, 31.29, 31.17, 31.09, 31.17, 31.34,
     +  31.61, 31.87, 32.62, 33.16, 33.84, 34.62, 35.85, 36.67, 37.51,
     +  38.72, 39.52, 40.29, 41.00, 41.97, 42.54, 43.01, 43.38, 43.82,
     +  44.01, 44.14, 44.24, 44.34, 44.38, 44.40, 44.40, 44.36, 44.33,
     +  44.28, 44.21, 44.07, 43.96, 43.84, 43.71, 43.48, 43.33, 43.17,
     +  43.02, 42.79, 42.64, 42.49, 42.36, 42.16, 42.05, 41.94, 41.85,
     +  41.71, 41.63, 41.54, 41.46, 41.35, 41.27, 41.20, 41.13, 41.02,
     +  40.95, 40.88, 40.81, 40.72, 40.66, 40.61, 40.57, 40.49, 40.44,
     +  40.37, 40.30, 40.19, 40.11, 40.04, 39.92, 39.84, 39.76, 39.69,
     +  39.57, 39.51, 39.45, 39.40, 39.32, 39.26, 39.21, 39.16, 39.09,
     +  39.05, 39.01, 38.98, 38.93, 38.90, 38.88, 38.85, 38.82, 38.81,
     +  38.79, 38.77, 38.75, 38.73, 38.72, 38.70, 38.68, 38.66, 38.65,
     +  38.64, 38.62, 38.60, 38.59, 38.57, 38.55, 38.53, 38.51, 38.49/
       data sigpb /
     + 505.81,485.41,470.03,455.70,441.22,420.60,407.94,395.19,377.36,
     + 366.03,355.16,345.07,329.94,320.60,311.79,303.19,291.17,283.77,
     + 276.68,270.09,260.15,254.24,248.43,243.19,234.87,229.91,225.10,
     + 220.26,213.46,209.12,204.98,200.98,195.36,191.56,187.98,184.49,
     + 179.41,176.04,172.86,169.81,165.19,162.28,159.35,156.55,152.41,
     + 149.95,147.59,145.11,141.66,139.59,137.45,135.41,132.43,130.53,
     + 128.59,126.81,124.18,122.54,120.90,119.34,117.13,115.67,114.22,
     + 112.04,110.54,109.13,107.70,105.63,104.23,102.97,101.64, 99.62,
     +  98.30, 96.99, 95.81, 93.97, 92.88, 91.80, 90.69, 89.03, 87.97,
     +  86.92, 85.89, 84.36, 83.38, 82.39, 81.40, 79.89, 78.92, 77.95,
     +  76.96, 75.50, 74.59, 73.68, 72.81, 71.51, 70.73, 69.95, 69.19,
     +  68.07, 67.39, 66.71, 66.04, 65.13, 64.54, 63.99, 63.42, 62.55,
     +  62.05, 61.52, 61.02, 60.23, 59.74, 59.27, 58.79, 58.07, 57.64,
     +  57.20, 56.77, 56.10, 55.66, 55.23, 54.41, 53.90, 53.42, 52.94,
     +  52.20, 51.73, 51.30, 50.85, 50.20, 49.79, 49.47, 49.16, 48.72,
     +  48.46, 48.15, 47.84, 47.42, 47.17, 46.92, 46.68, 46.31, 46.12,
     +  45.90, 45.70, 45.42, 45.28, 45.13, 44.99, 44.80, 44.68, 44.58,
     +  44.48, 44.35, 44.28, 44.21, 44.15, 44.06, 44.00, 43.95, 43.90/
       data sigkp /
     +   9.05,  8.86,  8.71,  8.58,  8.45,  8.28,  8.19,  8.11,  8.02,
     +   7.97,  7.94,  7.92,  7.90,  7.91,  7.92,  7.95,  8.01,  8.05,
     +   8.10,  8.17,  8.28,  8.37,  8.47,  8.56,  8.74,  8.86,  8.98,
     +   9.13,  9.35,  9.51,  9.66,  9.83, 10.10, 10.27, 10.45, 10.63,
     +  10.89, 11.07, 11.24, 11.40, 11.65, 11.81, 11.96, 12.12, 12.37,
     +  12.54, 12.71, 12.87, 13.19, 13.43, 13.70, 13.99, 14.49, 14.87,
     +  15.28, 15.71, 16.37, 16.79, 17.19, 17.57, 18.05, 18.30, 18.49,
     +  18.68, 18.74, 18.76, 18.75, 18.69, 18.63, 18.57, 18.50, 18.39,
     +  18.34, 18.29, 18.25, 18.18, 18.14, 18.09, 18.05, 17.98, 17.92,
     +  17.87, 17.81, 17.72, 17.66, 17.62, 17.59, 17.55, 17.52, 17.49,
     +  17.46, 17.44, 17.43, 17.42, 17.41, 17.41, 17.40, 17.40, 17.40,
     +  17.40, 17.40, 17.40, 17.40, 17.41, 17.41, 17.42, 17.42, 17.43,
     +  17.43, 17.44, 17.44, 17.45, 17.45, 17.46, 17.47, 17.48, 17.48,
     +  17.48, 17.49, 17.49, 17.50, 17.50, 17.50, 17.49, 17.49, 17.49,
     +  17.49, 17.50, 17.50, 17.50, 17.51, 17.52, 17.54, 17.55, 17.58,
     +  17.60, 17.62, 17.63, 17.66, 17.68, 17.70, 17.73, 17.76, 17.79,
     +  17.81, 17.83, 17.87, 17.90, 17.92, 17.95, 17.99, 18.01, 18.03,
     +  18.06, 18.10, 18.13, 18.16, 18.19, 18.23, 18.26, 18.29, 18.33/
       data sigkm /
     +  76.56, 74.23, 72.44, 70.86, 69.28, 66.99, 65.66, 64.35, 62.47,
     +  61.33, 60.23, 59.18, 57.65, 56.68, 55.73, 54.78, 53.21, 52.17,
     +  51.31, 50.48, 49.27, 48.55, 47.89, 47.30, 46.48, 45.91, 45.31,
     +  44.74, 43.78, 43.13, 42.48, 41.77, 40.57, 39.85, 39.01, 38.22,
     +  37.05, 36.34, 35.66, 35.01, 34.23, 33.98, 33.73, 33.53, 33.26,
     +  33.24, 33.41, 33.68, 34.13, 34.53, 34.94, 35.52, 36.37, 37.02,
     +  37.58, 38.00, 38.11, 37.88, 37.44, 36.88, 35.77, 35.12, 34.34,
     +  33.26, 32.60, 31.91, 31.29, 30.54, 30.13, 29.77, 29.45, 29.03,
     +  28.80, 28.58, 28.35, 27.97, 27.70, 27.40, 27.10, 26.70, 26.48,
     +  26.29, 26.07, 25.73, 25.52, 25.31, 25.13, 24.86, 24.72, 24.57,
     +  24.42, 24.17, 24.01, 23.89, 23.77, 23.60, 23.49, 23.38, 23.29,
     +  23.14, 23.05, 22.97, 22.90, 22.79, 22.71, 22.64, 22.56, 22.45,
     +  22.37, 22.28, 22.20, 22.07, 21.99, 21.91, 21.82, 21.69, 21.61,
     +  21.53, 21.45, 21.34, 21.27, 21.20, 21.11, 21.06, 21.00, 20.95,
     +  20.87, 20.83, 20.80, 20.77, 20.71, 20.68, 20.65, 20.62, 20.57,
     +  20.54, 20.52, 20.50, 20.46, 20.43, 20.40, 20.38, 20.35, 20.32,
     +  20.30, 20.28, 20.25, 20.24, 20.23, 20.21, 20.20, 20.20, 20.19,
     +  20.19, 20.18, 20.18, 20.18, 20.17, 20.16, 20.15, 20.14, 20.12/
       data sigpi /
     +   5.02,  5.69,  6.20,  6.81,  7.56,  9.29, 10.54, 12.24, 15.11,
     +  17.54, 20.45, 23.24, 28.67, 32.75, 37.60, 42.78, 52.32, 59.55,
     +  67.49, 75.47, 87.11, 94.75,101.54,106.41,111.56,113.00,112.73,
     + 110.52,105.96,101.22, 95.06, 88.38, 77.77, 70.51, 63.86, 57.83,
     +  49.65, 45.30, 41.45, 38.23, 34.33, 32.40, 30.86, 29.68, 28.76,
     +  28.61, 28.53, 28.41, 28.04, 27.81, 27.89, 28.24, 29.65, 31.04,
     +  32.53, 33.91, 34.81, 34.92, 35.12, 35.32, 35.83, 36.30, 36.71,
     +  36.93, 36.79, 36.49, 36.18, 35.80, 35.59, 35.35, 35.12, 34.65,
     +  34.26, 33.85, 33.46, 32.97, 32.75, 32.61, 32.51, 32.37, 32.25,
     +  32.12, 31.96, 31.67, 31.47, 31.25, 31.02, 30.67, 30.44, 30.21,
     +  29.98, 29.65, 29.44, 29.23, 29.03, 28.72, 28.52, 28.33, 28.14,
     +  27.85, 27.66, 27.47, 27.30, 27.04, 26.90, 26.76, 26.65, 26.51,
     +  26.43, 26.35, 26.28, 26.19, 26.13, 26.07, 26.02, 25.92, 25.85,
     +  25.78, 25.70, 25.58, 25.51, 25.43, 25.32, 25.25, 25.18, 25.11,
     +  25.01, 24.94, 24.87, 24.81, 24.72, 24.65, 24.59, 24.54, 24.46,
     +  24.41, 24.36, 24.32, 24.25, 24.21, 24.18, 24.14, 24.09, 24.06,
     +  24.03, 24.00, 23.96, 23.93, 23.91, 23.88, 23.84, 23.82, 23.80,
     +  23.78, 23.76, 23.74, 23.73, 23.73, 23.72, 23.72, 23.72, 23.73/
c
c             25.49675=1./log(1.04)
c
      x      = 25.49675 * log(p * 10.)
      n      = x
      diff   = x - n
      n      = n + 1
      if (n .le. 0) n = 1
      if (n .gt. 162) n = 162
c
c             interpolate in sig tables.
c
      s      = sig(n,itype)
      if (n.eq.1 .or. n.eq.162) go to 50
      s = s + diff * (sig(n + 1,itype) - s)
  50  continue
c
c             now find abul-magd parameters
c
      is     = s
      if (is.le.0) is = 1
      if (is.ge.100) is = 99
      sb     = sigb(is) + (s - is) * (sigb(is + 1) - sigb(is))
      al     = alpha(is) + (s - is) * (alpha(is + 1) - alpha(is))
      ghsig  = sb*a**al
      return
      end
