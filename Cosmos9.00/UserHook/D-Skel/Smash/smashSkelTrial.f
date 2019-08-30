#if defined (KEKB) || defined (KEKA)
#define AVOIDFOOL
#endif
      implicit none
!   This version runs very quick with torelable accuracy : i.e,
!   each skeleton will get almost the same sum of energy and the number of particles
!  (example of energy sum list for 8000 cpu's is in Trial.rsl)
!
!
!           We smash the skeleton into NCPU skeeltons; but flesh (MCPU+MARGIN) skeletons
!           and modify the  fleshed one by multiplying a factor NCPU/MCPU
!           to the total number of particles. 
!           We select (MCPU+MARGIN) skeletons among NCPU skeletons randomly.
!
!
!       read skelton data and store each children as a complete
!       track so that each can be put to stack area as incidnet
!       particle.
!    smashed skeleton data format
!  Assume ncpu cpu's; For each cpu,  smashed skeleton files will be
!
!      skeleton001
!        cumnum, num, ir, Zfirst
!        Np
!        observed ptcles 1
!        observed ptcles 2
!
!        observed ptcles Np
!        nlowp  
!        track-1
!        track-2
!        ... 
!        track-nlowp
!     other skeleton file( skeleton002,...)
!        cumnum, num, ir, Zfirst
!        0
!        nlowp  
!        track-1
!        track-2
!        ...
!        track-nlowp
!        ....
!
#include "Ztrack.h"
#include "Zearth.h"
      include "../../SkelFlesh/Zprivate.h"
      include "Zprivate2.h"


      type(child):: cc
      integer icon
      integer klena
      integer i, nlow, cumnum, num, ir(2)
      integer ii
      type(track):: Zfirst
      character*120  outdir
      character*120 skelfilebase
      character*100  basename
      character*100  filename
      character*100  input
      character*100  hostlist, thinhostlist
      character*15  field(3)
      integer n, j, k, nr, ll, kgetenv2

      hostlist = ' '
      thinhostlist = ' '
      
      ll = kgetenv2("NCPU", msg)
      read( msg(1:ll), *) Ncpu
      ll = kgetenv2("MCPU", msg)
      read( msg(1:ll), *) Mcpu
      ll = kgetenv2("MARGIN", msg)
      read( msg(1:ll), *) Margin


      ll = kgetenv2("SKELETON", msg)
      skelfilebase=msg(1:ll)
      ll = kgetenv2("SMSKELDIR", msg)
      outdir = msg(1:ll)
      ll = kgetenv2("SKELNAME", msg)
      basename= msg(1:ll)
      ll = kgetenv2("HOSTLIST", msg)
      if(ll .gt. 0) hostlist = msg(1:ll)

      ll = kgetenv2("THINHOSTLIST", msg)
      if(ll .gt. 0) thinhostlist = msg(1:ll)

!        binary open
      call copenfw2(11, skelfilebase, 2, icon)
      if(icon .ne. 1) then
         write(msg,*) skelfilebase(1:klena(skelfilebase)),
     *    ' could not be opened properly'
         call cerrorMsg(msg, 0)
      endif
      write(msg,*) "# of cpu's=",Ncpu, Mcpu
      call cerrorMsg(msg, 1)
      if(Ncpu .lt. 1 .or. Ncpu .gt. MaxCPU) then
         call cerrorMsg("# of cpu's > MaxCPU <1 ",0)
      endif

!        open  output smashed skeleton files
      k = klena(outdir)
      if(  outdir(k:k) .ne. '/') then
         k = k + 1
         outdir(k:k)= '/'
      endif
      write(msg, '(a,a)') 'output directory is ',
     *   outdir(1:k)
      call cerrorMsg(msg, 1)
      write(msg,*) Ncpu,
     *    ' files will be created there as '//
     *    basename(1:klena(basename))//'0001 etc' 
      call cerrorMsg(msg, 1)

!
      
      if(hostlist .ne. ' ') then
         call copenf(12, hostlist, icon)
         if(icon .ne. 0 ) then
            call cerrorMsg(hostlist, 1)
            call cerrorMsg(' could not be opened', 0)
         endif
      else
         write(0,*) ' hostlist not given'
         stop 1234
      endif

      if( Mcpu .le. Ncpu ) then
         if( thinhostlist .ne. ' ' ) then
            call copenf(13, thinhostlist, icon)
            if(icon .ne. 0 ) then
               call cerrorMsg(thinhostlist, 1)
               call cerrorMsg(' could not be opened', 0)
            endif
         else
            write(0,*) ' thinhostlist not given'
            stop 1234
         endif
      else
         write(0,*) ' Mcpu =', Mcpu, ' > Ncpu=', Ncpu
         stop 1234
      endif

      do i = 1, Ncpu
         read(12, '(a)') input
!            input may be like:      1  hosta   2.5
         field(1) = ' '
         field(2) = ' '
         field(3) = ' '
         call ksplit(input, 30, 3, field,  nr)
         read(field(1), '(i6)' )  numba(i)
         if(  numba(i)  .gt. Ncpu .or. numba(i) .lt. 1 ) then
            write(0,*) i,  '-th line has ivalid number=',
     *           numba(i), ' in Hosts'
            stop 0000
         endif
         if(nr .le. 2) then
            cpupw(i) = 1.0
         else
            read(field(3), * )  cpupw(i)
         endif
      enddo
      close(12)

      if(Mcpu .le. Ncpu) then
         do i = 1, Mcpu+Margin
            read(13, '(a)') input
!            input may be like:      1  hosta   2.5
            field(1) = ' '
            field(2) = ' '
            field(3) = ' '
            call ksplit(input, 30, 3, field,  nr)
            read(field(1), '(i6)' )  ii
            if(  ii  .gt. Ncpu .or. ii .lt. 1 ) then
               write(0,*) i,  '-th line has ivalid number=',
     *            ii, ' in ThinHosts'
               stop 0000
            endif
            numba(ii) = -ii
         enddo
         close(13)
      endif

      do i = 1, Ncpu
         if( numba(i) .lt. 0 ) then
           write(filename,'(a,i5.5)') 
     *        basename(1:klena(basename)), i
           skelfile(i) = outdir(1:klena(outdir))//filename
!            We don't open files here; too  many files
!            might not  be opened simultaneously.
!         call copenfw2(basefn+i, skelefile, 2, icon)
!         if(icon .ne. 0) then
!            call cerrorMsg(skelfile(i), 1)
!            call cerrorMsg('could not be opened properly',1)
!            call cerrorMsg('maybe they already exist', 0)  
         endif
      enddo


!      ------------


      do while(.true.)
         read(11, end=100) cumnum, num, ir, 
#if defined (AVOIDFOOL)
#include "ZavoidUnionMap.h"
#else
     *   Zfirst
#endif
         do i = 1, Ncpu
            if( numba(i) .lt. 0 ) then
!                 open one file at a time and
!                 if possible.
#if defined (AVOIDFOOL)
               open(basefn+i, file=skelfile(i), position="append",
#else
               open(basefn+i, file=skelfile(i), access="append" , 
#endif
     *              form="unformatted",   iostat=icon)
               if(icon .eq. 0) then
                 write(basefn+i)  cumnum, num, ir, 
#if defined (AVOIDFOOL)
#include "ZavoidUnionMap.h"
#else
     *            Zfirst
#endif

               else
                  write(0,*)  skelfile(i), " could not be opened"
                  stop 0000
               endif
               close(basefn+i)
            endif
         enddo

         read(11) Np
         call cerrorMsg('------------', 1)
         write(msg, *) Np, ' ptcls are observed ones in skeleton'
         call cerrorMsg(msg, 1)
         if(Np .gt. Maxob) then
            call cerrorMsg(
     *      'It is too large; enlarge Maxob', 0)
         endif

         do i = 1, Np
            read(11) oo(i)
         enddo
         nlow = 1
         ctc=0
         do while (nlow .ge. 0)
            read(11) nlow, pp
!               nlow = 0, if pp.asflag=-1.
            do i = 1, nlow
               read(11) cc
               if(ctc .lt. Maxp) then
                  ctc = ctc + 1
                  call movetrack(cc, ct(ctc) )
               else
                  call cerrorMsg(
     *                 'too many particles in skeleton',1)
                  call cerrorMsg(
     *            'Enlarge Maxp in Zprivate2.h', 0)
               endif
            enddo
         enddo

         write(msg,*)
     *   '# of total ptcls at flesh=',ctc
         call cerrorMsg(msg, 1)

!             1 event data is ready now in oo and ct.
!             distribute particles to ncpu
!                  first sort ct by energy
!///////////
         write(0,*) ' sort by energy starts'
         call sortbyerg
         write(0,*) ' sort by energy ended'
!/////////////////
!             deploy particles to Ncpu so that
!             sum energy on each cpu is roughly  the same
         if(ctc .lt. Ncpu) then
            n = ctc
            write(msg, *) '# of ptcls < Ncpu'
            call cerrorMsg(msg, 1)
         else
            n = Ncpu
         endif
         write(0,*) ' distribute particles to', n, ' cpus'
         call distribute( n )

         write(0,*) ' starts to write sub skeletons'
         call memoforcpu( n )

         call issuemsg(  Ncpu )
      enddo

 100  continue
      call cerrorMsg('all events have been smashed',1)
!      do i = 1, Ncpu
!         close(basefn+i)
!      enddo

      end
!     ----------------------------
      subroutine distribute( n )
      implicit none
#include "Ztrack.h"
      include "../../SkelFlesh/Zprivate.h"
      include "Zprivate2.h"
      integer i, k
      integer n, j
      logical bigtosmall

      
      do i = 1, Ncpu
         sumergi(i)= 0.
         sumergw(i) = 0.
         nOnCpu(i) = 0
      enddo
      do i = 1, n 
!          max energy ptcl for i-th cpu
         sumergi(i) = erg(idx(i))
         sumergw(i) = erg(idx(i)) / cpupw(i)
         nOnCpu(i) = 1
         idxlist(1, i) = idx(i)
         idxlocal(i) = i
      enddo
!          if all cpupw =1, next two not needed
      call kqsortd(sumergw, idxlocal, n)
      call ksortinv(idxlocal, n)

!///////////
!      write(0,*) ' top E=',(sumergi(i), i=1, n)
!      write(0,*) ' idx=',(idx(i), i=1, n)
!////////
!          next explanation is for cpupw = 1
!            erg      idx     sumergi   idxlocal  nOnCpu  idxlist
!                                                          1,1
!        1    9         5        30       1        1       5
!        2    1         3        18       2        1       3
!    n   3   18         7        15       3        1       7
!        4    5         8
!        5   30         1
!        6    4         4
!        7   15
!        8   13      
!        .   
!        .   
!        .              6
!      ctc   .          2
!
! after j= 4
!    sumergi idxlocal nOnCpu  idxlist
!                              1   2
!     30      1        1       5
!     18      2        1       3   
!     28      3        2       7   8
!  after j=5
!    sumergi idxlocal nOnCpu  idxlist
!                              1   2
!     30      1        1       5
!     27      3        2       3   1
!     28      2        2       7   8
!  after j=6
!    sumergi idxlocal nOnCpu  idxlist
!                              1   2  3
!     30      1        1       5
!     32      3        3       3   1  4
!     28      2        2       7   8
!
      k = n + 1
      bigtosmall=.true.


      do j = n+1, ctc
         if( bigtosmall ) then
            k  = k - 1
            if( k .eq. 0 ) then
               bigtosmall = .false.
               k = 1
            endif
         else
            k = k + 1
            if( k .gt. n ) then
               bigtosmall = .true.
               k =n
            endif
         endif
                                                                                         
         nOnCpu( k )  =  nOnCpu( k )   + 1
         if( nOnCpu( k ) .gt.  MaxPtclPerCpu ) then
            write(msg, *)
     *      '# of ptcls on a cpu', k, ' is ',  nOnCpu( k ) ,
     *     ' exceeded limit=', MaxPtclPerCpu
            call cerrorMsg(msg, 1)
            call cerrorMsg('Enlarge MaxPtclPerCpu in Zprivate2.h',0)
         endif
         idxlist( nOnCpu(k), k ) = idx(j)
         sumergw(k) = sumergw(k) + erg(idx(j))/cpupw(k)
         sumergi(k) = sumergi(k) + erg(idx(j))
      enddo

      end
!     *************************
      subroutine memoforcpu( n ) 
      implicit none
#include "Ztrack.h"
      include "../../SkelFlesh/Zprivate.h"
      include "Zprivate2.h"
      
      integer n
      integer  navob, navobc, navobx
      integer i, j, icon
      real  avob, resob
      real*8  u
      integer cpuc
      type(track)::Zfirst  ! name is some reason   
                        
!          we distribute Np observed ptcls (at skeleton making time)
!          almost equally to Mcpu  cpu (not to Margin);
!          If some hosts fails,  observed ones will be lost since
!          hosts of Margin have no observed ptcls. This may not have
!          any dangour since, observed ptcls will be only at the core
!          region and we don't need such region for detector simulation
!        number of average ptcls 

      avob = Np
      avob = avob/Mcpu
!      navobx = max(Np/Ncpu, 1)
      navobx = Np/Mcpu
!      if( Np .eq. 0 ) navobx = 0

      resob = avob-navobx   !  0<= resob< 1

      navobc = 0
!/////////////////
      write(0,*) ' navobx=', navobx, ' resob=',resob
!///////////
      cpuc = 0
      do i = 1, Ncpu
         if(numba(i)  .lt. 0 ) then
            cpuc = cpuc + 1
            navob = navobx
            call rndc(u)
            if( u .lt.  resob) then
               navob = navobx + 1
            endif
            if(navobc+navob .gt. Np .or. cpuc .eq. Mcpu ) then
               navob = Np -  navobc
            endif
#if defined (AVOIDFOOL)
           open(basefn+i, file=skelfile(i), position="append",
#else
           open(basefn+i, file=skelfile(i), access="append",
#endif
     *      form="unformatted",      iostat=icon)
            if(icon .eq. 0) then
!/////////////
               write(0,*) ' cpu ',i, ' obs=',navob
!//////////////
              write(basefn+i) navob
!

               do j = navobc +1, navobc+navob
                  write(basefn+i) oo(j)
               enddo
               navobc = navobc + navob
            else
               write(0,*) ' skelfile=', skelfile(i), " cannot be opened"
               stop 11111
            endif

!   ***      enddo
!   ***      do i = 1, Ncpu
!cc         if(i .eq. 1) then
!              for the first skeleton, put observed ptcls
!c            write(basefn+i) Np
!c            do j = 1, Np
!c               write(basefn+1) oo(j)
!c            enddo
!c         else
!c            write(basefn+i) 0
!c         endif
            write(basefn+i)  nOnCpu(i)
            do j = 1, nOnCpu(i) 
#if defined (AVOIDFOOL)
              Zfirst =ct( idxlist(j, i)  )
              write(basefn+i) 
#include "ZavoidUnionMap.h"
#else
              write(basefn+i) ct( idxlist(j, i)  )
#endif
            enddo
            close(basefn+i)
        endif
      enddo
      end
      subroutine issuemsg( n ) 
      implicit none
#include "Ztrack.h"
      include "../../SkelFlesh/Zprivate.h"
      include "Zprivate2.h"
      
      integer n
      integer i

      msg = ' cpu#   cpuPW    Sum E        # of ptcls'
!      msg = 'cpu#     Sum E      # of ptcls'
      call cerrorMsg(msg, 1)
      do i = 1, n
         write(msg,'(i6, f7.1, g16.7, i9)')
!         write(msg,'(i3, g16.7, i9)')
     *      numba(i), cpupw(i), sumergi(i), nOnCpu(i)
!     *      i,  sumergi(i), nOnCpu(i)
         call cerrorMsg(msg, 1)
      enddo
      end

      
      subroutine sortbyerg
      implicit none
#include "Ztrack.h"
      include "../../SkelFlesh/Zprivate.h"
      include "Zprivate2.h"

      integer i

      averg = 0.
      do i = 1, ctc
         erg(i) = ct(i).p.fm.p(4)
         averg = averg + erg(i) 
      enddo
      call kqsortd(erg, idx, ctc)
!       high to low
      call ksortinv(idx, ctc)
      if(ctc .gt. 0.) then
!            average total energy on 1 cpu
         averg = averg/ctc  * ncpu
      else
         call cerrorMsg('no child',1)
         return
      endif
      if( erg(idx(ctc) ) .gt. averg*1.1 ) then
!          max energy is too large. issue
!          warning
         write(msg,*) 'WARGNING: max E=', erg(idx(i)),
     *   ' is > average total energy for 1 cpu=',
     *   averg
         call cerrorMsg(msg, 1)
      endif
      end

      

      subroutine movetrack(f, t)
      implicit none
#include "Ztrack.h"
#include "Zearth.h"
      include "../../SkelFlesh/Zprivate.h"
      include "Zprivate2.h"

      type(child):: f
      type(track):: t

      t.p.code = f.code
      t.p.subcode = f.subcode
      t.p.charge = f.charge
      t.p.fm.p(1) = f.fm(1)
      t.p.fm.p(2) = f.fm(2)
      t.p.fm.p(3) = f.fm(3)
      t.p.fm.p(4) = f.fm(4)
      t.p.mass  = f.mass
      t.pos.xyz.r(1)  = pp.posx
      t.pos.xyz.r(2)  = pp.posy
      t.pos.xyz.r(3)  = pp.posz
      
      t.pos.depth = pp.depth 
      t.pos.height = pp.height
      t.pos.colheight = pp.colHeight 
      t.t = pp.atime 
      t.where = pp.where
      t.pos.radiallen =
     *     Eradius + pp.height
      t.pos.xyz.sys = 'xyz'
      t.vec.w.sys = 'xyz'
      t.wgt = 1.0
      t.asflag = 0
!      t.user = pp.user
      call cresetDirec( t )
      end


