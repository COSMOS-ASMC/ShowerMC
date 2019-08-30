#include "csoftenPiK.f"
#include "cutptcl.f"
      subroutine xBgRun
      implicit none
#include  "Zmaxdef.h"
#include "Zglobalc.h"
#include "Zmanagerp.h"
#include "Ztrack.h"
#include "Ztrackp.h"
#include "Ztrackv.h"
#include "Zcode.h"
#include "Zheavyp.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
#include  "Zstackv.h"
#include "Zprivate.h"
#include "Zprivate2.h"
#include "../Zabsorb.h"
      integer id  ! input.  1 ==> aTrack is going out from
!                                 outer boundery.
!                           2 ==> reached at an observation level
!                           3 ==> reached at inner boundery.
      type(track):: aTrack

!      type(track):: inci
!////////////
      type(coord):: pdir, cdir
!////////////
      type(coord):: tetafai
      
      character*128 input
      character*64 dirstr
      real sr, dr, tempr
      integer i, j, k, m,  icon
      integer ansites
      save ansites
      integer iij, code
      integer i1, i2, ic
      integer ir,  ifai, l, ridx, faiidx
      real*8  E0, cosz
      real*8  fai0, fai, sint
      real*8  delta  
      integer reducedTime
      integer NN
      integer klena
      integer w2hl(MaxNoOfSites)
      real*8 r, Eloss, rinmu, cosang
      real*8 dedt, dedtF, rho, dist, disto, BinFai
      real*8 aa
      real*8 wx, wy, wz, temp
      real   za
      real  de, Ek, f, molu
      real  dt, tmin 
      real*8  cvh2den
      data BinFai/30./
      integer ldep
!     integer ndummy
      character*9 ptcln(4)
      data ptcln/"Photons", "Electrons","Muons", "hadron"/
      character*9 ptcl2(3)
      data ptcl2/"Electrons", "Muons","All"/
      real power(4)
      integer nstr
      data power/1.,1.,1., 1./
      real  power2(3)
      data power2/1.,1.,1./
      character*128 title
      character*96 evid(nsites)
      save evid
      real*8 cog, cog2, sumne,  obstimes, Savederg(5)
      real*8 firstcdepth, dd
      logical dosort
      real*4  wt, stime
      real*8 sumEbydEdx, sumEbyDeath,sumEbyDeathNeu,sumEbyDeathNut
      real*8 sumEbyDeathE, sumEbyDeathG, sumEbyDeathMuPiK, 
     *      sumEbyDeathP, sumEbyDeathO
      real*8 sumEcrash, sumEspace
      real*8 sumAll, sumdEinAir, sumMissing, sumUncertain
      integer vn/2/ ! version number for the fnoB output
      save 
!/////////////
      real*8 pabs, rcore, sina, cs, sn,  cf, mom(3), Ek8, u
!/////////////

!     ***********************
      include "interface1.h"
!     *********************


      do i = 1, nsites
         w2hl(i) = 0
      enddo

      do i = 1, nsites
!             histdep(i) is the layer number
         if(histdep(i) .eq. 0)  exit
         ansites = i
         w2hl( histdep(i) ) = i
      enddo

      r=rmin
      dr = 10.**bin 

      do i = 1, nrbin
!            center of the bin:   
         rbin(i) = r
         r = r* dr
      enddo
#if defined (MACOSX)
#else
!            specify bin or ascii output
      call kwhistso( binw )
#endif

      return      
!     *********************************** hook for Beginning of  1 event
!     *  All system-level initialization for 1 event generation has been
!     *  eneded at this moment.
!     *  After this is executed, event generation starts.
!     *
      entry xBgEvent

      call cqIncident(inci, angle)
      E0 = inci.p.fm.p(4)
      if(inci.p.code .eq. kmuon) then
         call csetMuonPol(1.0d0)
      endif
      cosz = -angle.r(3)
      fai0 = atan2(-angle.r(2), -angle.r(1))*Todeg
      sint = sqrt(1.0-cosz**2)

      if(inci.p.code .eq. 9) then
         NN= inci.p.subcode
      elseif(inci.p.code .eq. 1) then
         NN=0
      else
         NN=1
      endif

      write(0,'("i ",  i6,  i4, g13.4,3f11.7,f7.1)')
     *   EventNo, inci.p.code,
     *   inci.p.fm.e,  -angle.r(1),  -angle.r(2), -angle.r(3)
       write(0,'(a, 1p, 6g15.5)')
     *    '### ', DetXaxis.r(1:3), DetZaxis.r(1:3)
      do i = 1, NoOfSites
         SumEloss(i) = 0.
         do j = 1, 4
            Ng(i) = 0.
            Ne(i) = 0.
            Nmu(i) = 0.
            Nhad(i) = 0.
         enddo
         dECent(i) = 0.
         do ifai = 1, nfai
            do ir= 1, nrbin
               dErfai(ir, ifai, i) = 0.
!               do j = 1, 4
!                  nrfaiAll(ir, ifai, j, i) = 0.
!               enddo
            enddo
         enddo
      enddo


!          estimate time minimum and time bin for eeach web sector
      
      do i = 1, ansites
         ldep = histdep(i)
         call cminTime2WebSec(ObsSites(ldep).pos.xyz,
     *        ldep, i,  webmin )
      enddo

#if defined (MACOSX)
#else
!     histogram: instanciate
!           t spectrum at each web sector
      if(tkrtspec) then
         do i = 1, ansites
            do j = 1, 4
!                at  center
               call kwhisti( tspec0(j, i),
     *          -5., 0.05, 200, b'00000')
               call kwhistai( tspec0(j, i),
     *         "Arrival time dist. of "//ptcln(j)//" at center",
     *         "t", "ptcls", .false., 0., "time", "ns")
!                   clear
               call kwhistc(tspec0(j, i))

               do ir=1, nrbin
                  do ifai=1, nfai
                     if(reducedTime .eq. 1) then
                        tmin = webmin(ir, 7, i)
                     else
                        tmin = webmin(ir, ifai,i)
                     endif
                     dt = 0.01*10.0**(bin*(ir-1))*100. ! approx core distnace m
                     dt = dt**0.75*1.e9/3.0e8/100. ! if sqrt 1m-->0.03 ns 10 m-->0.15 ns
                                       ! 100m 1ns 1km 5ns   4km 10ns
                                       ! dt**0.65  makes larger bin at large distance (<=x2)
                     if(j .eq. 4) dt=dt*10.0*ir/35.0 ! for delayed hadrons
                     dt= max(dt, 0.02)

                     call kwhisti( tspec(j, ir, ifai, i),
     *                    tmin,  dt, 2000,   b'00000')

                     call kwhistai( tspec(j, ir, ifai,  i),
     *                "Arrival time of "//ptcln(j)//" at (r,fai)",
     *                "rt", "ptcls", .false., 0., "time", "ns")
!                     clear 
                     call kwhistc(tspec(j, ir, ifai, i))
                  enddo
               enddo
            enddo
         enddo
      endif

!            lateral in each fai bin
      if(tkarspec) then

         do i = 1, ansites
            do j = 1, 4         ! g,e,mu,h
               do ifai = 1, nfai
                  call kwhisti(rspec(j, ifai, i),
     *                 rmin, bin, nrbin,  b'00011' )
                  call kwhistai(rspec(j, ifai, i), 
     *            "Lateral Dist. of "//ptcln(j)//" at  diff. azimuth",
     *            "ar", "ptcls", .true.,  power(j),   "r", "m.u")
!                     clear
                  call  kwhistc( rspec(j, ifai, i) )
               enddo
            enddo
         enddo
      endif
#endif

      obstimes = 0.

      return
!     ***************
      entry xObs(aTrack, id)
!
!     For id =2, you need not output the z value, because it is always
!     0 (within the computational accuracy).
!

      obstimes = obstimes + 1.d0
      if(mod(obstimes, 100000.d0) .eq. 0. ) then 
         dosort=.false.
         do i = 1, min(4,Stack_pos)
            if(Stack(i).p.fm.p(4) .ne. Savederg(i)) then
               Savederg(i)=Stack(i).p.fm.p(4) 
               dosort=.true.
            endif
         enddo
         if(dosort) then
            call csortStack
         endif
         write(0, *) ' obstimes=', obstimes, ' ptclE=',aTrack.p.fm.p(4)
         do i = 1, min(4,Stack_pos)
            write(0,*)' stack tops=', Stack(i).p.fm.p(4)
         enddo
      endif
!     ***************
      code = aTrack.p.code
      ldep =  aTrack.where
!     ************
      if(id .eq. 2 .and. code .le. 6 ) then  ! neglect rare ptcls
         wz = aTrack.vec.w.r(3) ! downgoing < 0
         if(wz .gt. 0) return
         wz = -wz
         r = sqrt( aTrack.pos.xyz.x**2 +
     *                 aTrack.pos.xyz.y**2 )
         molu = ObsSites(ldep).mu
         rinmu =r/molu
         sr = rinmu   ! single precision
         ridx = (log10( rinmu/rmin )/bin +0.5) +1 

         Ek = aTrack.p.fm.p(4) -aTrack.p.mass
         wt = aTrack.wgt  ! wt is single
         if(code .eq. kphoton) then
            Ng(ldep) = Ng(ldep) + aTrack.wgt
         elseif(code .eq. kelec) then
            Ne(ldep) = Ne(ldep) + aTrack.wgt
         elseif(code .eq. kmuon) then
            Nmu(ldep) = Nmu(ldep) + aTrack.wgt
         elseif(code .le. 6) then
            Nhad(ldep) = Nhad(ldep) +aTrack.wgt
         endif
!            ---------- compute energy loss 
         if(aTrack.p.charge .ne. 0  ) then
!             -----------------
!                      /|    |
!                     / |   1g/cm2 
!                    /A |    |
!            -------------------
!                  / ptcl direction  
!         get energy loss when aTrack goes some distance
!         of which vertical gramage is 1g/cm2.
!         Gramage the particle travel is 
!         1/cos where cos is the cos of angle (i.e, A if Fig)
!          in the detctor system.
!         1g/cm^2 = 10-3kg/10-4 m^2 =10 kg/m^2.
!         To travel  1 g/cm^2, the ptcl must
!         run dist kg/m^2
            if(abs(wz) .gt. 1.d-2) then
               dist =10./wz    ! in kg/m2/(g/cm2)
            else
!                 for safety
               dist =1000.
            endif

            call cqElossRate(dedt,dedtF) !  loss rate GeV/(kg/m^2)
!                     dedtF is the full eloss ; dedt is the restricted
!                     loss.  We may better use full here.
!                       energy in 1 g/cm2 of vertical direction
            Eloss =min( real(dedtF*dist), Ek)   !  GeV/(g/cm2)
            Eloss = Eloss*aTrack.wgt !  GeV/(g/cm2)
            SumEloss(ldep)=SumEloss(ldep) + Eloss
         else
            Eloss=0.
         endif

         if(code .ge. 4) code=4
         if( aTrack.p.charge .ne. 0   .or. 
     *       w2hl(ldep) .gt. 0 ) then 
!                  fai
!              fai is  in    -15 to 345  (for dfai=30.)                                     
            aa=atan2(aTrack.pos.xyz.y, aTrack.pos.xyz.x)*
     *        Todeg -fai0
            fai = aa/Todeg
            aa= mod(aa + 360.d0,   360.d0)
            if(aa .gt. (360.d0-dfai/2.0d0)) aa= aa-360.d0
            faiidx=(aa+dfai/2.0d0) /dfai + 1
            if(ridx .ge. 1 .and. ridx .le. nrbin) then
               dErfai(ridx, faiidx, ldep) = dErfai(ridx, faiidx, ldep)  
     *             +  Eloss         
            elseif(ridx .le. 0) then
               dECent(ldep) = dECent(ldep) + Eloss
            endif
! 
!  do following   for specified histo layers (typically only 1 layer)
!
#if defined  (MACOSX)
#else
            if( w2hl(ldep)  .gt. 0 ) then
               i = w2hl(ldep)
               if(tkarspec) then
                  call kwhist( rspec(code, faiidx,  i), 
     *            sr, wt)
               endif

               if( tkrtspec ) then
                  stime =  aTrack.t 
                  if(reducedTime .eq. 1) then
                     delta =  r*(cos(fai) + 1.)*sint*1.d9/c ! ns
                     stime = stime + delta
                  endif
                  ir = ridx
                  if(ir .lt. 1) then
                     call kwhist( tspec0(code, i), 
     *                  stime, wt)
                  elseif(ir .le. nrbin) then
                     call kwhist( tspec(code, ir, faiidx,  i), 
     *                    stime, wt)
                  endif
               endif
            endif
#endif
         endif
      endif
      return
!     **************
      entry xEnEvent
!     **************
!        replace  @ # % in basefilename by hostname, etc
!        and put it in basefilename2
!        
      
      write(0,*) 'ev#=',EventNo,
     *   ' generation phase finished. now writing data'

      call cgetfname(basefilename, basefilename2)
      call cqFirstID(firstcdepth)
      firstcdepth = firstcdepth* 0.1     ! in g/cm2  First col depth.



      if(ObserveAS) then
         cog = 0.
         sumne = 0.
         do i = 1, NoOfASSites
            if(i .gt. 1 .and. i  .lt. NoOfASSites ) then
               dd =(ASDepthList(i+1) - ASDepthList(i-1))/2.0
            elseif(i .eq. 1) then
               dd =(ASDepthList(2) - ASDepthList(1))
            else
               dd =(ASDepthList(NoOfASSites) -
     *              ASDepthList(NoOfASSites-1))
            endif
            cog = cog + ASObsSites(i).esize*dd*ASDepthList(i)
            sumne= sumne +ASObsSites(i).esize*dd
         enddo
!          0.1 is for g/cm2
         cog = cog*0.1/sumne

         cog2 = 0.
         sumne = 0.
         do i = 1, NoOfASSites
            if( ASObsSites(i).age .gt.
     *          (2.0-ASObsSites(NoOfASSites).age))  then
               if(i .gt. 1 .and. i  .lt. NoOfASSites ) then
                  dd =( ASDepthList(i+1) - ASDepthList(i-1))/2.0
               elseif(i .eq. 1) then
                  dd =(ASDepthList(2) - ASDepthList(1))
               else
                  dd =(ASDepthList(NoOfASSites) -
     *              ASDepthList(NoOfASSites-1))
               endif
               dd = dd
               cog2 = cog2 + ASObsSites(i).esize*ASDepthList(i)*dd
               sumne= sumne +ASObsSites(i).esize*dd
            endif
         enddo
         if(sumne .gt. 0.) then
            cog2 = cog2*0.1/sumne
         else
!              to deep penetration
            cog2 = ASDepthList(NoOfASSites)*0.1
         endif

         filename = basefilename2(1:klena(basefilename2))//".hyb"
         call copenfw2(fnoB, filename, 1, icon)

         write(fnoB,
     *   '("h ", i4,  3i3, 1pE11.3, 0p 3f11.7, 1pE11.3, 0p,
     *     2f7.0,i2,a )')
     *      EventNo,  inci.p.code,
     *      inci.p.subcode, inci.p.charge,
     *      inci.p.fm.e, -angle.r(1), -angle.r(2), -angle.r(3),
     *      firstcdepth, cog, cog2, vn, ' /'

         sumEbydEdx = 0.
         sumEbyDeathG =0.
         sumEbyDeathE =0.
         sumEbyDeathMuPiK =0.
         sumEbyDeathNeu = 0.
         sumEbyDeathNut = 0.
         sumEbyDeathP = 0.
         sumEbyDeathO = 0.
         sumEbyDeath = 0.
         sumUncertain = 0.
         sumEcrash = 0.
         sumEspace = 0.

         do i = 1, NoOfASSites 
            if(Eabsorb(1) .ne. 0) then
               write(fnoB, '("t ", i3, 2f7.1,  2f6.3,
     *         1p14g12.3 )')
     *           i, 
     *          ASDepthList(i)*0.1,  ASObsSites(i).mu,
     *          ASObsSites(i).age,   ASDepthList(i)*0.1/cog2, 
     *          Ng(i), Ne(i), Nmu(i), Nhad(i),
     *          ASObsSites(i).esize, SumEloss(i), 
     *          dEbydEdx(i), dEbyDeath(i), 
!                   next ones are from 7.51
     *          dEbyDeathG(i),  dEbyDeathE(i), dEbyDeathMuPiK(i), 
     *          dEbyDeathP(i),  dEbyDeathNut(i), dEbyDeathO(i)

               if(i .le. Eabsorb(2) ) then
!                    to see E consv. we should not count
!                    level > Eabsorb(2).
                  sumEbydEdx = sumEbydEdx + dEbydEdx(i)
                  sumEbyDeath = sumEbyDeath + dEbyDeath(i)
                  sumEbyDeathNeu = sumEbyDeathNeu +dEbyDeathNeu(i)
                  sumEbyDeathNut = sumEbyDeathNut +dEbyDeathNut(i)
                  sumEbyDeathG = sumEbyDeathG + dEbyDeathG(i)
                  sumEbyDeathE = sumEbyDeathE + dEbyDeathE(i)
                  sumEbyDeathMuPiK = sumEbyDeathMuPiK +
     *                   dEbyDeathMuPiK(i)
                  sumEbyDeathP = sumEbyDeathP +dEbyDeathP(i)
                  sumEbyDeathO = sumEbyDeathO +dEbyDeathO(i)
               endif
            else
               write(fnoB, '("t ", i3, 2f7.1,  2f6.3,
     *         1p6E11.3 )')
     *           i, 
     *          ASDepthList(i)*0.1,  ASObsSites(i).mu,
     *          ASObsSites(i).age,   ASDepthList(i)*0.1/cog2, 
     *          Ng(i), Ne(i), Nmu(i), Nhad(i),
     *          ASObsSites(i).esize, SumEloss(i)
            endif
         enddo
         if(Eabsorb(1) .ne. 0) then
            do i = 1, 7
               sumEcrash = sumEcrash + Ecrash(i)
               sumEspace = sumEspace + Espace(i)
            enddo
            write(fnoB,'("b ", 1p7E11.3)') (Espace(i), i=1,7)
            write(fnoB,'("b ", 1p7E11.3, i4)') (Ecrash(i), i=1,7),
     *      Eabsorb(2)
            write(fnoB,
     *       '("c ",1p7E11.3)' ) 
     *      MaxEbreak, MaxRelEbreak, SumEdiff, SumAbsEdiff,
     *      MaxEbreak(1)/inci.p.fm.p(4)

            sumMissing =  sumEcrash + sumEspace + sumEbyDeathNeu
            sumUncertain = sumEbyDeathNut
            sumdEinAir =  sumEbydEdx + sumEbyDeath
            sumAll = sumdEinAir + sumMissing + sumUncertain

            write(fnoB,'("s ", 1p8E11.3)') 
     *       sumEbydEdx, sumEbyDeath, sumdEinAir,
     *       sumEcrash, sumEspace, sumEbyDeathNut,
     *       sumEbyDeathNeu, sumAll

            write(fnoB,'("r ", 1p4E11.3)') 
     *      sumdEinAir/E0, sumUncertain/E0, sumMissing/E0, sumAll/E0
!                normalized one
            write(fnoB,'("n ", 1p4E11.3)') 
     *      sumdEinAir/sumAll, sumUncertain/sumAll,sumMissing/sumAll,
     *      1.0
!                additional info for more details
            write(fnoB,'("a ", 1p5g12.3 )')
     *      sumEbyDeathG,  sumEbyDeathE,  sumEbyDeathMuPiK,
     *      sumEbyDeathP,  sumEbyDeathO
         endif
         write(fnoB,*)
         close(fnoB)
      endif

      do i = 1, ansites
         j=histdep(i)
         write(evid(i), 
     *   '(i3, i5,  f5.2, f5.2,
     *   f7.1,  i4)')
     *   histdep(i), int( ASDepthList(j)*0.1 ),  
     *   ASObsSites(j).age, ASDepthList(j)*0.1/cog2,
     *   ASObsSites(j).mu, int(cog2)
      enddo
#if defined  (MACOSX)
#else
      if(tkarspec) then
         filename = basefilename2(1:klena(basefilename2))//"-r.hist"
         call copenfw2(fnoL, filename, binw, icon)
        do i = 1, ansites
            k=histdep(i)
            do j = 1, 4
               write(dirstr,'(a,"/d",i4, "/")')
     *              ptcln(j), int( DepthList(k)*0.1 )
               call kseblk(dirstr, "|", nstr)
               do l = 1, nfai
                  call kwhistdir(rspec(j, l,  i),  dirstr)
!                  call kwhists(  rspec(j, l, i), 0. )
                  call kwhists(  rspec(j, l, i), 0. )
                  call kwhistev( rspec(j, l, i), EventNo)
                  call kwhistid( rspec(j, l,  i), evid(i))
                  call kwhistp( rspec(j, l, i),  fnoL)
!                        *********** deallocate ********                            
                  call kwhistd( rspec(j, l, i) )
               enddo  ! code loop                                                   
            enddo ! fai loop                                                        
         enddo ! depth loop      

         close(fnoL)
      endif

      if( tkrtspec ) then
         filename = basefilename2(1:klena(basefilename2))//"-t.hist"
         call copenfw2(fnoT, filename, binw, icon)

         do i = 1, ansites
            do j = 1, 4
               call kwhists( tspec0(j,i), 0.)
               call kwhistev( tspec0(j,i), EventNo)
               call kwhistid( tspec0(j,i), evid(i))
               k=histdep(i)
               dirstr = " "
               write( dirstr,'(a,"/d",i4, "/")')
     *              ptcln(j), int( ASDepthList(k)*0.1 )
               call kseblk( dirstr, "|", nstr)
               call kwhistdir( tspec0(j,i),  dirstr )
               call kwhistp( tspec0(j,i),  fnoT )
!                 *********** deallocate ********         
               call kwhistd( tspec0(j,i) )
            enddo
         enddo

         do i = 1, ansites
            do j = 1, 4
               do ifai= 1, nfai
                  do ir= 1, nrbin
                     call kwhists( tspec(j,ir, ifai,i), 0.)
                     call kwhistev(tspec(j,ir, ifai,i), EventNo)
                     call kwhistid( tspec(j,ir, ifai,i), evid(i))
                     dirstr = " "
                     write(dirstr,'(a,"/d",i4, "/F",i2,"/")')
     *                    ptcln(j), int( DepthList(k)*0.1), ifai
                     call kseblk(dirstr, "|", nstr)
                     call kwhistdir(tspec(j,ir, ifai,i),  dirstr)
                     call kwhistp( tspec(j,ir, ifai,i),  fnoT)
!                        *********** deallocate ********                            
                     call kwhistd( tspec(j, ir, ifai, i) )
                  enddo
               enddo
            enddo  ! code loop                                                      
         enddo   ! depth loop    
         close(fnoT)
      endif
#endif
!          output web data
      if(tkweb) then
         filename = basefilename2(1:klena(basefilename2))//".nrfai"
         call copenfw2(fnoN, filename, 1, icon)
         
         write(fnoN,
     *   '(i4,1pE11.3, 0p,i3, f8.4, 1pE11.3,0p, 4i4, 1p,8g11.3)')
     *   EventNo, E0, NN, cosz, firstcdepth, nrbin, nfai, ansites,
     *   NoOfSites, KEminObs   ! this is not exist in the older version          
!                                                                      
!           dE/dx lateral                                              
         do i = 1, NoOfSites
            do k = 1, nfai
               write(fnoN, '("dE/dx",f7.1, 3i4)' )
     *            DepthList(i)*0.1, i, i, k
               write(fnoN, '(1p10E11.3)')
     *             ( dErfai(m,k,i), m=1,nrbin  ), dECent(i)
!                                same  center value is put for all fai
            enddo
         enddo
         close(fnoN)
      endif

      write(0,*) 'ev#=',EventNo,' finished completely'

      end
