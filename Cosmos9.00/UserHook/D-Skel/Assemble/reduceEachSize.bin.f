#define  FNODATDEF 33
#include "./realLimit.h"

!     reduce the file size; exmine the number of particles (N)
!     in each (r,fai) bin and if it exceeds a limit,
!     each partilce is discarded with the probabilty of
!     1-limit/N.
!       input file: combined  .nrfai data.  give it by environmental
!                   varialbe NRFAIFILE
!                   each  .dat file (binary file)    main output from sim. 
!                   give it by env. vari.  DATFILE
!                   env. vari.  HEADER  for the first set  give it "yes"
!                   else 'no'.  
!       output file:  size reduced file corresponding to .dat; this is put  stdout (ascii)
!                   so this must be concatinated to get combied .dat file.
!                   nrfaifile for each .dat; -r is attached as .nrfai-r.
!                    you must combine .nrfai-r files;
!                   rm old each .nrfai file. and rename .nrfai-r to .nrfai
!                   and use assemNrfai.csh in Assemble.
!  ********"all" data in   .nrfai-r  is really the "all" so this one should not be
!           added any more in  assemNrfai
#include "../FleshHist/asinfo.f"
#include "../FleshHist/asdensity.f"
#include "../FleshHist/crecprob.f"
      implicit none
#include "Zmaxdef.h"
#include "Zobs.h"
#include "../FleshHist/Zprivate0.h"
#include "../FleshHist/Zprivate1.h"
#include "../FleshHist/Zprivate3.h"

      integer  ldep,  
     *        charge, ridx, faiidx, codex

      real*8 u
      integer ndepth
      integer EventNo
      integer*2 code
      real*8   Et,  wx, wy, wz
      real    E0
      parameter (ndepth= nsites)
      real nrfaiRec0(nrbin, nfai, 4, ndepth)
      real rnrfaiRec(nrbin, nfai, 4, ndepth)
      real pnrfaiRec(nrbin, nfai, 4, ndepth)
      real nrfaiAll0(nrbin, nfai, 4, ndepth)
      real dErfai(nrbin, nfai,  ndepth)
      real dErfai2(nrbin, nfai,  ndepth)
      real recprob(nrbin, 4, ndepth)
      character*20 field(15)
      integer EvNo0
      integer l0(ndepth), recl0(ndepth)
      real limit(4), Emin
      real rat, all, rec
      integer klena
      real intdep(ndepth), recdep(ndepth)



      integer nrbina, nfaia, ansites0
      integer leng, i, j, k, l, leng2, lengh
      integer maxsites, nr
      integer i0, j0, k0, reci0, recj0, reck0
      integer NN
      integer icon0, iconx, icont
      character*128 input0
      character*100 nrfaifile, nrfaifile2, datfile, datfile2
      character*3 id
      character*5 id5
!      character*5 header
      integer fnonrfai, fnonrfai2,  fnodat2
      real  cosz, age, sum, Nx, depth, wgt, wwgt, prob
      real nptcls(nrbin, 4,  ndepth)
      integer indivdep(ndepth),  packeddepidx(ndepth), depidx
!      real rnptcls(nrbin, 4,  ndepth)
      integer  icon, kgetenv2, nrec
      logical SeeLowdE, accept
      integer newfmt
      fnonrfai=11
      fnonrfai2=21
!      fnodat = 6
!      fnodat2 = 8
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
      leng = kgetenv2("SeeLowdE", datfile)
      if(leng .le. 0)  then
         write(0,*) ' SeeLowdE not given'
         stop
      endif
      SeeLowdE = datfile(1:leng) .eq. "yes"

      leng2 = kgetenv2("DATFILE", datfile)
      if(leng2 .le. 0) then
         write(0,*) "Env. DATFILE not given"
         stop 11112
      endif
      call copenfw2(fnodat, datfile, 2, icon)
      if(icon .ne. 1) then
         write(0,*) ' input main data not exist'
         write(0,*) ' file=',datfile
         stop 1113
      endif

      datfile2 =" "
      datfile2 = datfile(1:leng2)//"-r"
      nrfaifile2 =" "
      nrfaifile2 = datfile(1:leng2-3)//"nrfai-r"
      call copenfw2(fnonrfai2, nrfaifile2, 1, icon)
      if(icon .ne. 0) then
         write(0,*) ' error ', nrfaifile2, ' cannot be created'
         stop 1235
      endif
!      lengh= kgetenv2("HEADER", header)

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
               Emin = 500.d-6
               newfmt = 0
            elseif(nr .eq. 9 ) then
               read(input0(1:klena(input0)), *) 
     *           EvNo0, E0,NN, cosz, limit(1), nrbina, nfaia, ansites0,
     *           maxsites 
               Emin = 500.d-6
               newfmt = 1
            elseif(nr .eq. 13) then
               read(input0(1:klena(input0)), *) 
     *           EvNo0, E0,NN, cosz, Emin, nrbina, nfaia, ansites0,
     *           maxsites, limit 
               newfmt = 2
            else
               write(0,*) ' header fmt strage '
               write(0,*) input0
               stop 01234
            endif

!            write(0,*) input0

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
     *               id, recdep(i), recl0(i), i0, j0, k0
                     indivdep(i)=recl0(i)
                     packeddepidx(recl0(i))=i  !  original dep index to packed indx
                     
!
!        when the above is written
!        l = indivdep(i)
!        write(fnonrfai, '("rec",f7.1, 4i4)' )
!     *     ASDepthList(l)*0.1, l, i, j, k
                                   

                     if(i0 .ne. i .or. j .ne. j0 .or. k .ne. k0) then
                        write(0,*) ' intdep, i0,j0,k0=',
     *                      recdep(i), i0, j0, k0, ' strange'
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
            do i = 1, maxsites
               do j = 1, 4
                  do k = 1, nfai
                     read(fnonrfai, '(a,f7.1, 4i4)' )
     *                id,  intdep(i), l0(i), i0, j0, k0
                     indivdep(i)=l0(i)
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
!              dErfai
            do i = 1, maxsites
               do k = 1, nfai
                  read(fnonrfai, '(a,f7.1, 3i4)' )
     *                 id5,  intdep(i), l0(i), i0,  k0
                  indivdep(i)=l0(i)
                  if(i0 .ne. i .or.  k .ne. k0) then
                     write(0,*) ' intdep, i0,k0=',
     *                      intdep(i), i0, k0, ' strange'
                     stop 98765
                  endif
                  if( id5 .ne. "dE/dx") then
                     write(0,*) 'id=',id5,' strange'
                     stop 9999
                  endif
                  read(fnonrfai, *)
     *                 ( dErfai(l,k,i), l=1,nrbin )
                  if(SeeLowdE) read(fnonrfai, *)
     *                 ( dErfai2(l,k,i), l=1,nrbin )
               enddo
            enddo
         else
!      
            write(0,*) ' all nrfai data has been read'
            limit(1) = REALLIMITg
            limit(2) = REALLIMITe
            limit(3) = REALLIMITmu
            limit(4) = REALLIMITh
            if(newfmt .eq.  0) then
               write(0, '(i2, 1pE11.3,i3, 0pf7.1, 1pE11.3,0p, 3i4)')
     *        EvNo0, E0,NN, cosz, limit(1),  nrbina, nfaia, ansites0
            elseif( newfmt .eq. 1 ) then
               write(0, '(i2, 1pE11.3,i3, 0pf7.1, 1pE11.3, 0p,4i4)')
     *        EvNo0, E0,NN, cosz, limit(1),  nrbina, nfaia, ansites0,
     *        maxsites
            else
              write(0, '(i2, 1pE11.3,i3, 0pf7.1,1pE11.3,0p,4i4,4f8.9)')
     *        EvNo0, E0,NN, cosz, Emin,  nrbina, nfaia, ansites0,
     *        maxsites, limit
           endif               
         endif
      enddo
 1000 continue
      input0= ' '
!************* reset limit
!      limit=min(4000.0, limit)
      limit(1) = REALLIMITg
      limit(2) = REALLIMITe
      limit(3) = REALLIMITmu
      limit(4) = REALLIMITh

      if(newfmt .ne. 2) then
         write(0,*) ' limit=',limit(1), ' is being used'
      else
         write(0,*) ' limit=',limit, ' ares being used'
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
     *              recdep(i), indivdep(i), i, j, k, ' prob' 
               write(0, '(1p10E11.3)')
     *               (pnrfaiRec(l,k,j,i), l=1, nrbin)
!c/////
            enddo
         enddo
      enddo
!  -------------
!            main input data; header Obsolete
!      read(fnodat)  EventNo, code,   Et,  wx, wy, wz
!      if( header(1:lengh) .eq. "yes") then
!         write(*,'("i ",  i3,  i4, g13.4,3f11.7)')
!     *        EventNo, code,   Et,  wx, wy, wz
!      endif

      nrec= 0
      do while(.true.)
         read(fnodat,end=100) bufc, (buf(i), i=1, bufc)
!/////////////
!         write(0,*) ' bufc=',bufc
!///////////
         do i = 1, bufc
            nrec= nrec+1
            ldep = buf(i).ldep
            depidx = packeddepidx(ldep)
            faiidx=  buf(i).faiidx
            ridx = buf(i).ridx
            codex = min(buf(i).code, 4)
            wgt = buf(i).wgt
            prob = pnrfaiRec(ridx, faiidx, codex, depidx) 
            if( prob .gt. 1.) then
               wwgt = wgt
               accept = .true.
            else
               prob = prob * wgt
               if(prob .gt. 1.) then
                  accept = .true.
                  wwgt = prob
               else
                  call rndc(u)
                  if(u .lt. prob) then
                     accept = .true.
                     wwgt= 1.
                  else
                     accept = .false.
                  endif
               endif
            endif
            if(accept) then
               rnrfaiRec(ridx, faiidx, codex, depidx)=
     *              rnrfaiRec(ridx, faiidx, codex, depidx) + wwgt
               write(*,
     *  '(6i3, 1pE11.3, 0p,f6.1,1p2E11.3,0p, 2f8.4,f10.6,1pE11.3)')
     *         buf(i).ldep,  buf(i).code,  buf(i).subcode,
     *         buf(i).charge, buf(i).ridx, buf(i).faiidx,
     *         buf(i).rinmu, buf(i).fai, buf(i).Ek,
     *         buf(i).t, buf(i).wx, buf(i).wy, buf(i).wz,
     *         wwgt
            endif
         enddo
      enddo
 100  continue

      write( fnonrfai2,
     *     '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3,0p,4i4,4f8.0)' )
     *     EvNo0, E0, NN, cosz, Emin, nrbin, nfai, ansites0,
     *     maxsites, limit

      do i = 1, ansites0
         do j = 1, 4
            do k = 1, nfai
               l = recl0(i)
               write(fnonrfai2, '("rec",f7.1, 4i4)' )
     *          recdep(i), l, i, j, k
               write(fnonrfai2, '(1p10E11.3)')
     *             ( rnrfaiRec(l,k,j,i), l=1,nrbin )
            enddo
         enddo
      enddo
      do i = 1, maxsites
         do j = 1, 4
            do k = 1, nfai
!               l = indivdep(i)
               l = l0(i)
               write(fnonrfai2, '("all",f7.1, 4i4)' )
     *          intdep(i), l, i, j, k
               write(fnonrfai2, '(1p10E11.3)')
     *             ( nrfaiAll0(l,k,j,i), l=1,nrbin )
            enddo
         enddo
      enddo
!       dErfai
      do i = 1, maxsites
         do k = 1, nfai
!            l = indivdep(i)
            l = l0(i)
            write(fnonrfai2, '("dE/dx",f7.1, 3i4)' )
     *           intdep(i), l, i, k
            write(fnonrfai2, '(1p10E11.3)')
     *           ( dErfai(l,k,i), l=1,nrbin )
            if(SeeLowdE) write(fnonrfai2, '(1p10E11.3)')
     *           ( dErfai2(l,k,i), l=1,nrbin )
         enddo
      enddo
      
      write(0,*) 'end of run'
      write(fnonrfai2, *) 
      stop
 500  continue
      write(0,*) ' input error at record =', nrec
      end
