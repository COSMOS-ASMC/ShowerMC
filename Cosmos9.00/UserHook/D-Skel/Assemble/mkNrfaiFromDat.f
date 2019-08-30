      implicit none
!   input: h0  .nrfai file with invalid "rec" but
!               correct "all" and "dE/dx" part
!               give the name by env. NRFAIFILE
!          f0  .dat file. combined main output file
!               stdin
!   output:h1  .nrfai with correct "rec", "all" and
!              "dE/dx"'
!               stdout
!   this creates "rec" part of .nrfai data from
!  .dat file
!  

!
#include  "Zmaxdef.h"
#include "Zobs.h"
#include "../FleshHist/Zprivate0.h"

      integer ndepth
      parameter (ndepth= nsites)
      real nrfaiRec(nrbin, nfai, 4, ndepth)
      real nrfaiAll(nrbin, nfai, 4, ndepth)
      real dErfai(nrbin, nfai, ndepth)

      integer EvNo,EvNo0
      integer fn0
      integer kgetenv2
      integer klena
      real intdep(ndepth)
      real E0, cosz, limit(4) 
      integer  newfmt, nr, maxsites
      real Emin
      integer NN
      integer nrbina, nfaia, ansites0
      integer leng, i, j, k, l
      integer i0, j0, k0
      integer l0(ndepth)
      integer icon
      character*128 nrfai
      character*128 input
      integer code, subcode, charge, ldep
      integer w2il(ndepth)
      character*20 field(15)
      character*128 input0
      integer ridx, faiidx
      real  rinmu, fai, Ek, time, wx, wy, wz
      fn0 = 11
      leng = kgetenv2("NRFAIFILE", nrfai)
      call copenfw2(fn0, nrfai, 1, icon)
      if(icon .ne. 1)  then
         write(0,*) nrfai(1:leng)
         if( icon .eq. 0) then
            write(0,*) 'not exists'
         else
            write(0,*) ' cannot be opened '
         endif
         write(0,*) ' icon=',icon
         stop 9999
      else
         write(0,*)  nrfai(1:leng), ' opened'
      endif

      do while(.true.)
         input = ' '
         read( fn0, '(a)',  end=1000 ) input
         if(input .ne. " ") then

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
               Emin=500.D-6
               newfmt = 1
            elseif( nr .eq. 13) then
               read(input0(1:klena(input0)), *)
     *           EvNo0, E0,NN, cosz, Emin, nrbina, nfaia, ansites0,
     *           maxsites, limit
               newfmt = 2
            endif
            if(nrbina .ne. nrbin .or. nfaia .ne. nfai) then
               write(0,*)' nrbina=',nrbina, 'or  nfaia=',nfaia,
     *              ' differ from the def. in this prog'
               stop 5555
            endif
!           ********
            do i = 1, ansites0
               do j = 1, 4
                  do k= 1, nfai
                     read(fn0, '(3x, f7.1, 4i4)' )
     *                    intdep(i), l0(i), i0, j0, k0
                     w2il( l0(i) ) = i
                     if(i0 .ne. i .or. j .ne. j0 .or. k .ne. k0) then
                        write(0,*) ' intdep, i0,j0,k0=',
     *                      intdep(i), i0, j0, k0, ' strange'
                        stop 8888
                     endif
                     read(fn0, *)
     *                    ( nrfaiRec(l,k,j,i), l=1,nrbin )
                     do l = 1, nrbin
                        nrfaiRec(l,k,j,i) = 0.
                     enddo
                  enddo
               enddo
            enddo
!    ************
            do i = 1, ansites0
               do j = 1, 4
                  do k = 1, nfai
                     read(fn0, '(3x,f7.1, 4i4)' )
     *                 intdep(i), l0(i), i0, j0, k0
                     if(i0 .ne. i .or. j .ne. j0 .or. k .ne. k0) then
                        write(0,*) ' intdep, i0,j0,k0=',
     *                      intdep(i), i0, j0, k0, ' strange'
                        stop 9876
                     endif
                     read(fn0, *)
     *                    ( nrfaiAll(l,k,j,i), l=1,nrbin )
                  enddo
               enddo
            enddo
            do i = 1, ansites0
               do k = 1, nfai
                  read(fn0, '(5x,f7.1, 3i4)' )
     *                 intdep(i), l0(i), i0,  k0
                  if(i0 .ne. i .or.  k .ne. k0) then
                     write(0,*) ' intdep, i0,k0=',
     *                    intdep(i), i0,  k0, ' strange'
                     stop 9875
                  endif
                  read(fn0, *)
     *                 ( dErfai(l,k,i), l=1,nrbin )
               enddo
            enddo
         endif
      enddo
 1000 continue

!        nrfai for   1 event read
!            read .dat and fill "rec"
      do while(.true.)
         read(*, *, end=100)
     *        ldep,  code,  subcode,
     *        charge, ridx, faiidx,
     *        rinmu, fai, Ek, time, wx, wy, wz
         if(code .gt. 4 ) code = 4 
         i =  w2il(ldep)
         nrfaiRec(ridx, faiidx, code, i) =
     *        nrfaiRec(ridx, faiidx, code, i) + 1.0
      enddo
 100  continue

!          outpute .nrfai
      if(newfmt .eq. 0 ) then
 0       write(*,
     *     '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3,0p,3i4)')
     *     EvNo, E0, NN, cosz, limit(1), 
     *     nrbina, nfaia, ansites0
      elseif( newfmt .eq. 1) then
 0       write(*,
     *     '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3,0p, 4i4)')
     *     EvNo, E0, NN, cosz, limit(1), 
     *     nrbina, nfaia, ansites0, maxsites
      elseif( newfmt .eq. 2) then
 0       write(*,
     *     '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3,0p, 4i4,4f10.0)')
     *     EvNo, E0, NN, cosz, Emin,
     *     nrbina, nfaia, ansites0, maxsites, limit
      endif
      do i = 1, ansites0
         do j = 1, 4
            do k = 1, nfaia
               write(*, '("rec",f7.1, 4i4)' )
     *              intdep(i), l0(i), i, j, k
               write(*, '(1p10E11.3)')
     *              ( nrfaiRec(l,k,j,i), l=1,nrbin )
            enddo
         enddo
      enddo
      do i = 1, ansites0
         do j = 1, 4
            do k = 1, nfaia
               write(*, '("all",f7.1, 4i4)' )
     *              intdep(i), l0(i), i, j, k
               write(*, '(1p10E11.3)')
     *              ( nrfaiAll(l,k,j,i), l=1,nrbin )
            enddo
         enddo
      enddo
    
      do i = 1, ansites0
         do k = 1, nfaia
            write(*, '("dE/dx",f7.1, 3i4)' )
     *           intdep(i), l0(i), i,  k
            write(*, '(1p10E11.3)')
     *           ( dErfai(l,k,i), l=1,nrbin )
         enddo
      enddo
      write(*, *) 

      write(0,*) 'all data in the event  processed '
      end
