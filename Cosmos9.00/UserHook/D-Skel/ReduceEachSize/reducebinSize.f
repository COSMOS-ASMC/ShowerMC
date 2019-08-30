#define  FNODATDEF 33
!     reduce the file size; exmine the number of particles (N)
!     in each (r,fai) bin and if it exceeds a limit,
!     each partilce is discarded with the probabilty of
!     1-limit/N.
!       input file: .nrfai data (assembled).  
!                    give it by environmental varialbe NRFAIFILE
!                   .dat file   ( main output from sim. )
!                    give it by environmental variable DATFILE
!       output file: size reduced ascii output file: stdout.
!                    and .nrfai file for that; name will be
!                    such that: let .dat file name is
!                    xxx.dat, then xxx.nrfai.

#include "../FleshHist/asinfo.f"
#include "../FleshHist/asdensity.f"
#include "../FleshHist/crecprob.f"
      implicit none
#include "../FleshHist/Zprivate0.h"
#include "../FleshHist/Zprivate1.h"
#include "../FleshHist/Zprivate3.h"

      type(buffer):: abuf

      integer  ldep,  code, subcode,
     *        charge, ridx, faiidx, codex
      real  rinmu, fai,  Ek, time,
     *              wx, wy, wz
      real*8 u
      integer ndepth
      real    E0
      parameter (ndepth= nsites)
      real nrfaiRec0(nrbin, nfai, 4, ndepth)
      real rnrfaiRec(nrbin, nfai, 4, ndepth)
      real pnrfaiRec(nrbin, nfai, 4, ndepth)
      real nrfaiAll0(nrbin, nfai, 4, ndepth)
      real recprob(nrbin, 4, ndepth)

      integer EvNo0

      real limit
      real rat, all, rec, prob
      integer klena
      real intdep(ndepth)
      integer nrbina, nfaia, ansites0
      integer leng, i, j, k, l
      integer i0, j0, k0, l0
      integer NN
      integer icon0, iconx, icont
      character*128 input0
      character*100 nrfaifile, nrfaifile2, datfile
      character*3 id
      integer fnonrfai, fnonrfai2
      real  cosz, age, sum, Nx, depth
      real nptcls(nrbin, 4,  ndepth)
      integer indivdep(ndepth),  packeddepidx(ndepth), depidx
!      real rnptcls(nrbin, 4,  ndepth)
      integer  icon, kgetenv2, nrec


      fnonrfai=11
      fnonrfai2=21
      leng = kgetenv2("NRFAIFILE", nrfaifile)
      if(leng .le. 0) then
         write(0,*) "Env. NRFAIFILE not given"
         stop 11111
      endif
      call copenfw2(fnonrfai, nrfaifile, 1, icon)
      if(icon .ne. 1) then
         write(0,*) ' error cannot open', nrfaifile
         stop 0000
      endif
      datfile=' '
      leng = kgetenv2("DATFILE", datfile)
      if(leng .le. 0) then
         write(0,*) "Env. DATFILE not given"
         stop 11112
      endif
      call copenfw2(fnodat, datfile, 2, icon)
      if(icon .ne. 1) then
         write(0,*) ' error cannot open', datfile
         stop 0002
      endif
      nrfaifile2 = datfile(1:leng-3)//"rnrfai"
      call copenfw2(fnonrfai2, nrfaifile2, 1, icon)
      write(0,*) ' reduced .nrfai is', trim(nrfaifile2)
      write(0,*) ' open icon=', icon
      do k = 1,  ndepth
         packeddepidx(k)=0
         do j = 1, 4
            do l = 1,nfai
               do i = 1, nrbin
                  rnrfaiRec(i, l, j, k)=0
               enddo
            enddo
         enddo
      enddo
!
      do while(.true.)
         input0 = ' '
         read(fnonrfai, '(a)',  end=1000 ) input0
         if(input0 .ne. " ") then
            read(input0(1:klena(input0)), *) 
     *           EvNo0, E0,NN, cosz, limit,  nrbina, nfaia, ansites0

            write(0,*) input0

            if(nrbina .ne. nrbin .or. nfaia .ne. nfai) then
               write(0,*)' nrbina=',nrbina, 'or  nfaia=',nfaia,
     *              ' differ from the def. in this prog'
               stop 5555
            endif
            if(ansites0 .gt. ndepth) then
               write(0,*) ' too many depths'
               stop 6666
            endif
!           ********
            do i = 1, ansites0
               do j = 1, 4
                  do k = 1, nfai
                     read(fnonrfai, '(a, f7.1,  4i4)' )
     *               id, intdep(i), l0, i0, j0, k0
                     indivdep(i)=l0
                     packeddepidx(l0)=i  !  original dep index to packed indx
                     
!
!        when the above is written
!        l = indivdep(i)
!        write(fnonrfai, '("rec",f7.1, 4i4)' )
!     *     ASDepthList(l)*0.1, l, i, j, k
                                   

                     if(i0 .ne. i .or. j .ne. j0 .or. k .ne. k0) then
                        write(0,*) ' intdep, i0,j0,k0=',
     *                      intdep(i), i0, j0, k0, ' strange'
                        stop 8888
                     endif
                     if( id .ne. "rec") then
                        write(0,*) 'id=',id, ' strange'
                        stop 5678
                     endif
!
                     read(fnonrfai, *)
     *                    ( nrfaiRec0(l,k,j,i), l=1,nrbin )
                  enddo
               enddo
            enddo
!    ************
            do i = 1, ansites0
               do j = 1, 4
                  do k = 1, nfai
                     read(fnonrfai, '(a,f7.1, 4i4)' )
     *                id,  intdep(i), l0, i0, j0, k0
                     indivdep(i)=l0
                     if(i0 .ne. i .or. j .ne. j0 .or. k .ne. k0) then
                        write(0,*) ' intdep, i0,j0,k0=',
     *                      intdep(i), i0, j0, k0, ' strange'
                        stop 9876
                     endif
                     if( id .ne. "all") then
                        write(0,*) 'id=',id,' strange'
                        stop 9999
                     endif
                     read(fnonrfai, *)
     *                    ( nrfaiAll0(l,k,j,i), l=1,nrbin )
                  enddo
               enddo
            enddo
         else
!      
            write(0,*) ' all nrfai data has been read'
            write(0, '(i2, 1pE11.3,i3, 0pf7.1, 1pE11.3, 3i4)')
     *       EvNo0, E0,NN, cosz, limit,  nrbina, nfaia, ansites0
         endif
      enddo
 1000 continue
      input0= ' '
!************* reset limit
!      limit=min(4000.0, limit)
      write(0,*) ' limit=',limit, ' is being used'
!************ 
      do i = 1,  ansites0
         do j = 1, 4
            do k = 1,nfai
               do l = 1, nrbin
                  if(nrfaiRec0(l, k, j, i) .gt. limit) then
!                       accept with this prob.
                     pnrfaiRec(l, k, j, i)=
     *                    limit/nrfaiRec0(l, k, j, i) 
                  else
                     pnrfaiRec(l, k, j, i)=1.0
                  endif
               enddo
!///////
!               write(0, '(f7.1,  4i4,a)' )
!     *              intdep(i), indivdep(i), i, j, k, ' prob' 
!               write(0, '(1p10E11.3)')
!     *               (pnrfaiRec(l,k,j,i), l=1, nrbin)
!c/////
            enddo
         enddo
      enddo
!  -------------
!            main input data; header
      call binreadHead
      nrec= 0

      do while(.true.)
         call binread1(abuf, icon)
         if(icon .ne.  0) exit
         nrec= nrec+1
         depidx = packeddepidx(abuf.ldep)

         if(depidx .le. 0) then
            write(0,*) ' should not happen. depidx=',depidx
            write(0,*) 'nrec=',nrec, 'icon=',icon, ' ldep=',
     *         abuf.ldep, ' code=',abuf.code
            stop 9875
         endif

         call rndc(u)
         codex=min(abuf.code, 4)
         if(u .lt.
     *    pnrfaiRec(abuf.ridx, abuf.faiidx, codex, depidx) )  then
            rnrfaiRec(abuf.ridx, abuf.faiidx, codex, depidx)=
     *        rnrfaiRec(abuf.ridx, abuf.faiidx, codex, depidx) + 1

          if(abuf.wz .lt. 0.999) then
            write(*,
     *           '(6i3, 1pE11.3, 0p,f6.1,1p2E11.3,0p, 2f8.4,f10.6)')
     *           abuf.ldep,  abuf.code, abuf.subcode,
     *           abuf.charge, abuf.ridx, abuf.faiidx,
     *           abuf.rinmu, abuf.fai,
     *           abuf.Ek, abuf.t,
     *           abuf.wx, abuf.wy, abuf.wz
           else
             write(*,'(6i3,1pE11.3, 0p, f6.1, 1p2E11.3, a)')
     *           abuf.ldep,  abuf.code, abuf.subcode,
     *           abuf.charge, abuf.ridx, abuf.faiidx,
     *           abuf.rinmu, abuf.fai,
     *           abuf.Ek, abuf.t,
     *           ' 0 0 1'
           endif
        endif
      enddo
 100  continue
      write( fnonrfai2,
     *     '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3,3i4)' )
     *     EvNo0, E0, NN, cosz, limit, nrbin, nfai, ansites0

      do i = 1, ansites0
         do j = 1, 4
            do k = 1, nfai
               l = indivdep(i)
               write(fnonrfai2, '("rec",f7.1, 4i4)' )
     *          intdep(i), l, i, j, k
               write(fnonrfai2, '(1p10E11.3)')
     *             ( rnrfaiRec(l,k,j,i), l=1,nrbin )
            enddo
         enddo
      enddo
      do i = 1, ansites0
         do j = 1, 4
            do k = 1, nfai
               l = indivdep(i)
               write(fnonrfai2, '("all",f7.1, 4i4)' )
     *          intdep(i), l, i, j, k
               write(fnonrfai2, '(1p10E11.3)')
     *             ( nrfaiAll0(l,k,j,i), l=1,nrbin )
            enddo
         enddo
      enddo
      write(0,*) 'end of run'
      end
!       ************************
      subroutine binreadHead
!
#include "../FleshHist/Zprivate1.h"
#include "../FleshHist/Zprivate3.h"
      real*8  Et, wx, wy, wz
      integer EventNo
      integer*2  code
      integer i, icon
      
      type(buffer):: abuf
      save

      read(fnodat)
     *        EventNo, code,   Et,  wx, wy, wz
      write(*,'("i ",  i3,  i4, g13.4,3f11.7)')
     *        EventNo, code,   Et,  wx, wy, wz
      
      bufc = 0
      return

!     ******************
      entry  binread1(abuf, icon)
!     ***************
      if(bufc .eq. 0) then
         read(fnodat,end=100, err=500) bufc, (buf(i), i=1, bufc)
!        ****************************
!          bufc usage was wrong for the first 10^19 eV event made by
!          olixxx.   it is 0~99999; while correct one should be
!                    1~100000.  we neglect data at 0 and 100000.
         if(bufc .eq. 100000) bufc=bufc-1
      endif
      if(bufc .gt. 0) then
         abuf= buf(bufc)
         bufc = bufc -1
         icon = 0
      else
         write(0, *) ' strange bufc'
         stop 99999
      endif
      return
 100  continue
      icon = 1
      return
 500  continue
      write(0,*) ' input data read err'
      stop
      end
