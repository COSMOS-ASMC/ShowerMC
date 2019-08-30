!    suppose you have .nrfai file with invalid "rec" component while
!    "all" and "dE/dx" part are correct.
!    Also you have size reduced .dat file.
!    This program read .dat file and count data in each (r, fai) bin
!    and make correct  .nrfai data.
!   The program structure is the same as reduceSize.f but we accept
!   all particles in the .dat
!       input file: .nrfai data.  give it by environmental
!                   varialbe NRFAIFILE
!                    .dat file    main output from sim.
!                       standard input.
!       output file:  rfaifile with .nrfai-r extension.        

#include "../FleshHist/asinfo.f"
#include "../FleshHist/asdensity.f"
#include "../FleshHist/crecprob.f"

      implicit none
#include "Zmaxdef.h"
#include "Zobs.h"
#include "../FleshHist/Zprivate0.h"
#include "./realLimit.h"
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
      real dErfai(nrbin, nfai, ndepth)
      real recprob(nrbin, 4, ndepth)

      integer EvNo0

      real limit(4), Emin
      integer maxsites, nr, newfmt
      real rat, all, rec, prob
      integer klena
      real intdep(ndepth)
      integer nrbina, nfaia, ansites0
      integer leng, i, j, k, l
      integer i0, j0, k0, l0
      integer NN
      integer icon0, iconx, icont
      character*128 input0
      character*100 nrfaifile, nrfaifile2
      character*3 id
      character*5 id5
      character*20 field(15)
      integer fnonrfai, fnonrfai2
      real  cosz, age, sum, Nx, depth
      real nptcls(nrbin, 4,  ndepth)
      integer indivdep(ndepth),  packeddepidx(ndepth), depidx
!      real rnptcls(nrbin, 4,  ndepth)
      integer  icon, kgetenv2, nrec


      fnonrfai=11
      fnonrfai2=21
      nrfaifile=" "
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
      nrfaifile2 =" "
      nrfaifile2 = nrfaifile(1:leng)//"-r"
      call copenfw2(fnonrfai2, nrfaifile2, 1, icon)
      if(icon .ne. 0) then
         write(0,*) ' error ', nrfaifile2, ' cannot be created'
         stop 1235
      endif

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
            call  ksplit(input0,  20, 15, field,  nr)
            if(nr .eq. 8) then
               read(input0(1:klena(input0)), *)
     *           EvNo0, E0,NN, cosz, limit(1), nrbina, nfaia, ansites0
               maxsites = ansites0
               Emin=500.d-6
               newfmt = 0
            elseif(nr .eq. 9 ) then
               read(input0(1:klena(input0)), *)
     *           EvNo0, E0,NN, cosz, limit(1), nrbina, nfaia, ansites0,
     *           maxsites
               Emin=500.d-6
               newfmt = 1
            elseif( nr .eq. 13) then
               read(input0(1:klena(input0)), *)
     *           EvNo0, E0,NN, cosz, Emin, nrbina, nfaia, ansites0,
     *           maxsites, limit
               newfmt = 2
            endif
            write(0,*) input0
            limit(1) = REALLIMITg
            limit(2) = REALLIMITe
            limit(3) = REALLIMITmu
            limit(4) = REALLIMITh

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
!             dErfai
            do i = 1, ansites0
               do k = 1, nfai
                  read(fnonrfai, '(a,f7.1, 3i4)' )
     *                 id5,  intdep(i), l0, i0, k0
                  indivdep(i)=l0
                  if(i0 .ne. i .or. k .ne. k0) then
                     write(0,*) ' intdep, i0,k0=',
     *                    intdep(i), i0, k0, ' strange'
                     stop 98765
                  endif
                  if( id5 .ne. "dE/dx") then
                     write(0,*) 'id=',id5,' strange'
                     stop 9999
                  endif
                  read(fnonrfai, *)
     *                    ( dErfai(l,k,i), l=1,nrbin )
               enddo
            enddo
         else
!      
            write(0,*) ' all nrfai data has been read'
            if(newfmt .eq. 0) then
               write(0, '(i2, 1pE11.3,i3, 0pf7.1, 1pE11.3,0p, 3i4)')
     *       EvNo0, E0,NN, cosz, limit(1),  nrbina, nfaia, ansites0
            elseif(newfmt .eq. 1) then
               write(0, '(i2, 1pE11.3,i3, 0pf7.1, 1pE11.3,0p, 4i4)')
     *       EvNo0, E0,NN, cosz, limit(1),  nrbina, nfaia, ansites0,
     *       maxsites
            else
               write(0, 
     7        '(i2, 1pE11.3,i3, 0pf7.1, 1pE11.3, 0p,4i4, 4f10.0)')
     *       EvNo0, E0,NN, cosz, Emin,  nrbina, nfaia, ansites0,
     *       maxsites, limit
            endif
         endif

      enddo
 1000 continue
      input0= ' '
!************* reset limit
!      limit=min(4000.0, limit)
      if(newfmt .ne. 2) then
         write(0,*) ' limit=',limit(1), ' is being used'
      else
         write(0,*) ' limit=',limit, ' are being used'
      endif
!************ 
      do i = 1,  ansites0
         do j = 1, 4
            do k = 1,nfai
               do l = 1, nrbin
                  if(nrfaiRec0(l, k, j, i) .gt. limit(j)) then
!                       accept with this prob.
                     pnrfaiRec(l, k, j, i)=
     *                    limit(j)/nrfaiRec0(l, k, j, i) 
                  else
                     pnrfaiRec(l, k, j, i)=1.0
                  endif
               enddo
!///////
               write(0, '(f7.1,  4i4,a)' )
     *              intdep(i), indivdep(i), i, j, k, ' prob' 
               write(0, '(1p10E11.3)')
     *               (pnrfaiRec(l,k,j,i), l=1, nrbin)
!c/////
            enddo
         enddo
      enddo
!  -------------
!            main input data; header; obsolete
!      read(*,'(a)') input0
!      write(*,'(a)') input0(1:klena(input0))
      nrec= 0
      do while(.true.)
         read(*,'(a)', end=100, Err=500) input0
         if( index(input0(1:2), "i") .gt. 0 ) cycle
!         read(*, *, end=100, Err=500)
         read(input0,*)
     *      ldep,  code, subcode,
     *      charge, ridx, faiidx,
     *      rinmu, fai,
     *      Ek, time,
     *      wx, wy, wz
         nrec= nrec+1
         depidx = packeddepidx(ldep)

         if(depidx .le. 0) then
            write(0,*) ' should not happen. depidx=',depidx
            write(0,*) ' ldep=',ldep, ' code=',code, 'nrec=',nrec
            stop 9875
         endif
         codex=min(code, 4)
!         if(u .lt.
!     *    pnrfaiRec(ridx, faiidx, codex, depidx) )  then
          rnrfaiRec(ridx, faiidx, codex, depidx)=
     *        rnrfaiRec(ridx, faiidx, codex, depidx) + 1
      enddo
 100  continue

      if( newfmt .eq. 0) then
         write( fnonrfai2,
     *     '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3,0p,3i4)' )
     *     EvNo0, E0, NN, cosz, limit(1), nrbin, nfai, ansites0
      elseif( newfmt .eq. 1) then
         write( fnonrfai2,
     *     '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3,0p,4i4)' )
     *     EvNo0, E0, NN, cosz, limit(1), nrbin, nfai, ansites0,
     *     maxsites 
      else
         write( fnonrfai2,
     *     '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3,0p, 4i4,4f10.0)' )
     *     EvNo0, E0, NN, cosz, Emin, nrbin, nfai, ansites0,
     *     maxsites, limit 
      endif
         

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
!     dErfai
      do i = 1, ansites0
         do k = 1, nfai
            l = indivdep(i)
            write(fnonrfai2, '("dE/dx",f7.1, 3i4)' )
     *           intdep(i), l, i, k
            write(fnonrfai2, '(1p10E11.3)')
     *           ( dErfai(l,k,i), l=1,nrbin )
         enddo
      enddo
      
      write(0,*) 'end of run'
      write(fnonrfai2,*)
      stop
 500  continue
      write(0,*) ' input error at record =', nrec
      read(*,'(a)') input0
      write(0,*) input0
      end
