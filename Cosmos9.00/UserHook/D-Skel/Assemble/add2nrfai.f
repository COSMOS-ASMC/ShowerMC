      implicit none
!   Suppose  we have n nrfai data
!     h1, h2, h3, ....hn
!   1) cp h1 h0
!   2) h0 + h2 --> h; mv  h h0
!   3) h0 + h3 --> h; mv  h h0
!   ..
!   n) h0 + hn --> h; mv  h h0
! 
!   This program add two nrfai data; h0 + hx--> h
!            environmental variable
!   file h0: NRFAIFILE0
!   file hx: NRFAIFILEX
!   file h:  NRFAIFILE
!   "all" is should not touched:  GETALL 
!         
!
#include "Zobs.h"
#include "../FleshHist/Zprivate0.h"

      integer ndepth
      parameter (ndepth= nsites)
      real nrfaiRec0(nrbin, nfai, 4, ndepth)
      real nrfaiAll0(nrbin, nfai, 4, ndepth)

      real nrfaiRecx(nrbin, nfai, 4, ndepth)
      real nrfaiAllx(nrbin, nfai, 4, ndepth)

      real nrfaiRect(nrbin, nfai, 4, ndepth)
      real nrfaiAllt(nrbin, nfai, 4, ndepth)
   
      real dErfai0(nrbin, nfai, ndepth)
      real dErfai02(nrbin, nfai, ndepth)
      real dErfaix(nrbin, nfai, ndepth)
      real dErfaix2(nrbin, nfai, ndepth)
      real dErfait(nrbin, nfai, ndepth)
      real dErfait2(nrbin, nfai, ndepth)

      integer EvNo0, EvNox
      integer fn0, fnx, fnt
      integer kgetenv2
      integer klena
      real intdep(ndepth), recdep(ndepth)
      real E0, cosz, limit(4)
      integer NN
      integer nrbina, nfaia, ansites0, ansitesx
      integer maxsites
      integer newfmt
      integer leng, i, j, k, l, nr
      integer i0, j0, k0, reci0, recj0, reck0
      integer l0(ndepth), recl0(ndepth)
      integer icon0, iconx, icont
      character*128 nrfai0, nrfaix, nrfait, flesher
      character*128 input0, inputx, inputt
      character*20 field(15)
      character*5 getall
      real*4 Emin
      logical SeeLowdE

      fn0 = 2
      fnx = 3
      fnt = 4
      leng = kgetenv2("GETALL", getall)
      getall=getall(1:leng)

      write(0,*) ' getall =', getall

      leng = kgetenv2("NRFAIFILE0", nrfai0)
      call copenfw2(fn0, nrfai0, 1, icon0)
      if(icon0 .ne. 1)  then
         write(0,*) nrfai0(1:leng)
         if( icon0 .eq. 0) then
            write(0,*) 'not exists'
         else
            write(0,*) ' cannot be opened '
         endif
         write(0,*) ' icon=',icon0
         stop 9999
      else
         write(0,*)  nrfai0(1:leng), ' opened'
      endif
      leng = kgetenv2("NRFAIFILEX", nrfaix)
      call copenfw2(fnx, nrfaix, 1, iconx)
      if(iconx .ne. 1)  then
         write(0,*) nrfaix(1:leng)
         if( iconx .eq. 0) then
            write(0,*) 'not exists'
         else
            write(0,*) ' cannot be opened '
         endif
         write(0,*) ' icon=',iconx
         stop 9999
      else
         write(0,*)  nrfaix(1:leng), ' opened'
      endif

      
      leng = kgetenv2("NRFAIFILET", nrfait)
      call copenfw2(fnt, nrfait, 1, icont)
      if(icont .ne. 0)  then
         write(0,*) nrfait(1:leng)
         write(0,*) ' cannot be opened '
         write(0,*) ' icon=',icont
         stop 9999
      else
         write(0,*)  nrfait(1:leng), ' opened'
      endif
      leng = kgetenv2("FLESHDIR", flesher)

      leng = kgetenv2("SeeLowdE", input0)
      if(leng .le. 0) then
         write(0,*) "SeeLowdE not given"
         stop
      endif
      SeeLowdE = input0(1:leng) .eq. "yes"

         




      do while(.true.)
         input0 = ' '
         read( fn0, '(a)',  end=1000 ) input0
         inputx = ' '
         read( fnx, '(a)'  ) inputx 


         if(input0 .ne. " ") then
            call  ksplit(input0,  20, 15, field,  nr)
!//////////////
!            write(0,*) ' input0=', input0
!            write(0,*) ' nr=',nr
!///////////
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
            else
               write(0,*) ' header fmt strage '
               write(0,*) input0
               stop 01234
            endif
            if(nrbina .ne. nrbin .or. nfaia .ne. nfai) then
               write(0,*)' nrbina=',nrbina, 'or  nfaia=',nfaia,
     *              ' differ from the def. in this prog'
               stop 5555
            endif
            if(newfmt .eq. 0) then
               read(inputx(1:klena(inputx)), *) 
     *           EvNox, E0,NN, cosz, limit(1),
     *            nrbina, nfaia, ansitesx
            elseif(newfmt .eq.  1) then
               read(inputx(1:klena(inputx)), *) 
     *           EvNox, E0,NN, cosz, limit(1), nrbina, nfaia, 
     *           ansitesx, maxsites
            elseif(newfmt .eq. 2) then
               read(inputx(1:klena(inputx)), *) 
     *           EvNox, E0,NN, cosz, Emin, nrbina, nfaia, 
     *           ansitesx, maxsites, limit
            endif
               
            if(nrbina .ne. nrbin .or. nfaia .ne. nfai) then
               write(0,*)' nrbina=',nrbina, 'or  nfaia=',nfaia,
     *              ' differ from the def. in this prog'
               stop 6666
            endif
            if(ansites0 .ne. ansitesx) then
               write(0,*) ' ansites0=', ansites0,
     *              ' ansitesx=', ansitesx, ' diff '
               stop 7777
            endif

!           ********
            do i = 1, ansites0
               do j = 1, 4
                  do k= 1, nfai
                     read(fn0, '(3x, f7.1, 4i4)' )
     *                    recdep(i), recl0(i), reci0, recj0, reck0

                     if(reci0 .ne. i .or. j .ne. recj0 .or. 
     *                  k .ne. reck0) then
                        write(0,*) ' recdep, reci0,recj0,reck0=',
     *                  recdep(i), reci0, recj0, reck0, ' strange'
                        stop 8888
                     endif
                     read(fn0, *)
     *                    ( nrfaiRec0(l,k,j,i), l=1,nrbin )

                  enddo
               enddo
            enddo

!    ************
            do i = 1, maxsites
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
     *                    ( nrfaiAll0(l,k,j,i), l=1,nrbin )
                  enddo
               enddo
            enddo


            do i = 1, maxsites
               do k = 1, nfai
                  read(fn0, '(5x,f7.1, 3i4)' )
     *                 intdep(i), l0(i), i0,  k0
                  if(i0 .ne. i .or.  k .ne. k0) then
                     write(0,*) ' intdep, i0,k0=',
     *                    intdep(i), i0,  k0, ' strange'
                     stop 9875
                  endif
                  read(fn0, *)
     *                 ( dErfai0(l,k,i), l=1,nrbin )
                  if(SeeLowdE) read(fn0, *)
     *                 ( dErfai02(l,k,i), l=1,nrbin )
               enddo
            enddo


!           ************
!           ************
            do i = 1, ansites0
               do j = 1, 4
                  do k= 1, nfai
                     read(fnx, '(3x, f7.1, 4i4)' )
     *                    recdep(i), recl0(i), reci0, recj0, reck0
                     if(reci0 .ne. i .or. j .ne. recj0 .or.
     *                   k .ne. reck0) then
                        write(0,*) ' recdep, reci0,recj0,reck0=',
     *                      intdep(i), reci0, recj0, reck0, ' strange'
                        stop 9999
                     endif
                     read(fnx, *)
     *                    ( nrfaiRecx(l,k,j,i), l=1,nrbin )
                  enddo
               enddo
            enddo

!          ************            
            do i = 1, maxsites
               do j = 1, 4
                  do k = 1, nfaia
                     read(fnx, '(3x,f7.1, 4i4)' )
     *                 intdep(i), l0(i),  i0, j0, k0
                     if(i0 .ne. i .or. j .ne. j0 .or. k .ne. k0) then
                        write(0,*) ' intdep, i0,j0,k0=',
     *                      intdep(i), i0, j0, k0, ' strange'
                        stop 9877
                     endif
                     read(fnx, *)
     *                    ( nrfaiAllx(l,k,j,i), l=1,nrbin )
                  enddo
               enddo
            enddo

            do i = 1, maxsites
               do k = 1, nfaia
                  read(fnx, '(5x,f7.1, 3i4)' )
     *                 intdep(i), l0(i),  i0,  k0
                  if(i0 .ne. i .or.  k .ne. k0) then
                     write(0,*) ' intdep, i0,k0=',
     *                      intdep(i), i0, k0, ' strange'
                     stop 98778
                  endif
                  read(fnx, *)
     *                 ( dErfaix(l,k,i), l=1,nrbin )
                  if(SeeLowdE) read(fnx, *)
     *                 ( dErfaix2(l,k,i), l=1,nrbin )
               enddo
            enddo

!          ************            
         else

!            1 event read
            if(inputx .ne. " ") then
               write(0,*) ' event differ', inputx
               stop 1111
            endif
            do i = 1, ansites0
               do j = 1, 4
                  do k= 1, nfai
                     do l=1, nrbin
                        nrfaiRect(l,k,j,i) = 
     *                   nrfaiRec0(l,k,j,i) +
     *                   nrfaiRecx(l,k,j,i)
                     enddo
                  enddo
               enddo
            enddo
            if(getall .eq. 'yes') then
               do i = 1, maxsites
                  do j = 1, 4
                     do k= 1, nfai
                        do l=1, nrbin
                           nrfaiAllt(l,k,j,i) =
     *                          nrfaiAll0(l,k,j,i) +
     *                          nrfaiAllx(l,k,j,i)
                        enddo
                     enddo
                  enddo
               enddo
            endif    
            if(getall .eq. 'yes') then
               do i = 1, maxsites
                  do k= 1, nfai
                     do l=1, nrbin
                        dErfait(l,k,i) = 
     *                    dErfai0(l,k,i) +
     *                    dErfaix(l,k,i)
                        if(SeeLowdE) dErfait2(l,k,i) = 
     *                    dErfai02(l,k,i) +
     *                    dErfaix2(l,k,i)
                     
                     enddo
                  enddo
               enddo
            endif
            if( newfmt .eq. 0 ) then
               write(fnt,
     *       '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3, 0p, 3i4)')
     *        EvNox, E0, NN, cosz, limit(1), 
     *        nrbina, nfaia, ansites0
            elseif( newfmt .eq. 1) then
               write(fnt,
     *       '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3,0p, 4i4)')
     *        EvNox, E0, NN, cosz, limit(1), 
     *        nrbina, nfaia, ansites0, maxsites
            elseif( newfmt .eq. 2) then
               write(fnt,
     *       '(i2,1pE11.3, 0p,i3, f8.4, 1pE11.3,0p, 4i4,4f10.0)')
     *        EvNox, E0, NN, cosz, Emin,
     *        nrbina, nfaia, ansites0, maxsites, limit
            endif   

            do i = 1, ansites0
               do j = 1, 4
                  do k = 1, nfaia
                     write(fnt, '("rec",f7.1, 4i4)' )
     *                    recdep(i), recl0(i), i, j, k
                     write(fnt, '(1p10E11.3)')
     *                    ( nrfaiRect(l,k,j,i), l=1,nrbin )
                  enddo
               enddo
            enddo

            do i = 1, maxsites
               do j = 1, 4
                  do k = 1, nfaia
                     write(fnt, '("all",f7.1, 4i4)' )
     *                    intdep(i), l0(i), i, j, k
                     if(getall .eq. "yes") then
                        write(fnt, '(1p10E11.3)')
     *                    ( nrfaiAllt(l,k,j,i), l=1,nrbin )
                     else
                        write(fnt, '(1p10E11.3)')
     *                    ( nrfaiAll0(l,k,j,i), l=1,nrbin )
                     endif
                  enddo
               enddo
            enddo

            do i = 1, maxsites
               do k = 1, nfaia
                  write(fnt, '("dE/dx",f7.1, 3i4)' )
     *                 intdep(i), l0(i), i,  k
                  if(getall .eq. "yes") then
                     write(fnt, '(1p10E11.3)')
     *                    ( dErfait(l,k,i), l=1,nrbin )
                     if(SeeLowdE) write(fnt, '(1p10E11.3)')
     *                    ( dErfait2(l,k,i), l=1,nrbin )
                  else
                     write(fnt, '(1p10E11.3)')
     *                    ( dErfai0(l,k,i), l=1,nrbin )
                     if(SeeLowdE) write(fnt, '(1p10E11.3)')
     *                    ( dErfai02(l,k,i), l=1,nrbin )
                  endif
               enddo
            enddo

            write(fnt, *) 
         endif
      enddo
 1000 continue

      write(0,*) 'all data in the event  processed '
      end


            
