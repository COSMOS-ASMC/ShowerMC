#include "ZcosmosBD.h"
!               this is generic test routine for interaction code
!         inside Comsos.  h-A/ A'-A collision can be tested.
!        
          implicit none

#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zevhnp.h"
#include  "Zevhnv.h"
#include  "Zmass.h"
#include  "Zmanagerp.h"
#include  "Zmanager.h"
#include  "Ztrackp.h"
#include  "Zprivate.h"

          external cblkManager
          external cblkEvhnp
          integer code, subcode, charge, klena
          real*8 roots, pabs
          character*1  NULL

          integer jcode(maxDecay)

          type(ptcl):: pj, tg, pj2, a(nmax), cmspj, cmstg
          type(ptcl):: compj, comtg
          type(ptcl):: compjcms, comtgcms, cmspsave
!             projectile and target information (both befor
!             and after collision ) in different system.
!
          integer k, icon, ntp, j, nevent, nuccharge, ntp0
          real*8 x, y, eta, ek1, ek
          real*8 xcms,  ek1cms
          integer NEVTS,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IEMU
          real*8 EPN
          real*4  x2, y2, z2
          real*4  ext/10./
          type(ptcl):: labp

          real*8  energy
          real*8  gamma, beta
          real*8  crossint
          real pt, ptx, pty, ke, teta, pz, pzcms, ptcms
          real  kecms, tetacms

          character*100  msg
          character*100  tracefile
          character*80  prefix
          character*16 uid
          character*16 input
          integer ia, iz, leng, kgetenv2
          logical lroots, cms, xbyp
          integer dcy, i, kinc
!////////////
          real*8  tgtA, tgtZ, xs
!/////////////

!
          type(fmom):: rest
!         OmegaBaryon Dmeson  eta      gzai   lambda    lambdac  pi0  sigma
          data jcode/kbomega, kdmes, keta, kgzai,  klambda, klambdac,
     *              kpion, ksigma/

          call creadParam(5)
          call cmkSeed(InitRN(1), InitRN)
          call rnd1i(InitRN)     ! random number init.

          rest%p(1) = 0.
          rest%p(2) = 0.
          rest%p(3) = 0.
          rest%p(4) = 0.
!
!           get additional parameters written in user hookc
!
          call cqUHookc(1, msg)
          read(msg, *) lroots, code, subcode, charge, energy,
     *              ia, iz, cms, xbyp
          call cqUHookc(2, input)
          call cqUHookc(3, prefix)
          write(ErrorOut,*) input
!
!UserHooki = 0,       0,      0,       0,      0,        0,       0    0
!         OmegaBaryon Dmeson  eta      gzai   lambda    lambdac  pi0  sigma
!
          do i = 1, kindmx
             Jdecay(i) = 0
          enddo

          do i = 1, maxDecay
             dcy = 0
             call cqUHooki( i, dcy )
             Jdecay(jcode(i)) = dcy
          enddo

          if(DestEventNo(2) .eq. 0) then
             nevent =abs( DestEventNo(1))
          else
             nevent = abs(DestEventNo(2))
          endif

          if(Trace .eq. 0) then
             write(0,
     *      '(a,l2, a,i2, a,i3, a,i2, a,g13.3, a,i3, a,i2,'//
     *      ' a,l2, a,l2, a, a, a)')
     *      '# "roots=',lroots, ' proj%code=',code, ' subcode=',
     *      subcode, ' chg=',charge, '  E=',energy,' target A=',
     *      ia, ' Z=',iz, ' cms=',cms, ' xbyp=',xbyp,
     *      ' intmodel=', IntModel(1:klena(IntModel)),
     *      '"' 
            write(0, '(a)') 
     * '# "X" "Xcms" "Pz" "Pzcms" "y" "eta" "pt" "code" '//
     * ' "subcd" "chg" "mult" "teta" "tetacms" "K%E" "K%Ecms" "ev#"'
!           defalut trace.  fix the dirctor
         else
          if(TraceDir .eq. ' ') then
             call cgetLoginN(uid)
             TraceDir = '/tmp/'//uid(1:klena(uid))
          endif
         endif
!
!c             make incident
          call cmkptc(code, subcode,charge,  pj)
!            cms is for h-p
          if(lroots) then
             roots = energy
             if(pj%code .ne. kgnuc) then
                energy = ( roots**2 - pj%mass**2 - masp**2) /
     *                   (2*masp)
             else
                energy = ( roots**2/(2*masp) -  masp ) * subcode
             endif
          endif

!              set projectile energy and momentum
          pj%fm%p(1) =0.
          pj%fm%p(2) =0.
          pj%fm%p(4) = energy
          pj%fm%p(3) = sqrt(pj%fm%p(4)**2-pj%mass**2)
!
!           to form a CMS, once make proton target and h or p proj.
          call cmkptc(knuc, -1, 1,  tg)
          tg%fm = rest
          tg%fm%p(4) = tg%mass
          pj2 = pj
          if(pj%code .eq. kgnuc) then
             pj2%code = knuc
             pj2%charge = 1
             pj2%subcode = -1
             pj2%fm%p(4) = energy/subcode
             pj2%mass = masp
             pj2%fm%p(3) = sqrt(pj2%fm%p(4)**2-pj2%mass**2)
          endif

!             initialize. we are using as if from Comsos
          if(input .eq. ' ') then
             call cintModels('cosmos')
          else
             if(prefix .eq. ' ' ) then
                leng = kgetenv2("EPICSTOP", prefix)
                if(leng .eq. 0) then
                   call cerrorMsg('EPICSTOP is not given', 0)
                endif
!                        the last / is needed.
                prefix = prefix(1:leng)//'/Data/Media/'
             endif
             call cfixPrefix(prefix)
             call cintModels('check')
             if( index(IntModel, 'dpmjet3') .gt. 0) then
                call copenf(
     *           TempDev, input(1:klena(input))//".inp", icon)
                if(icon .ne. 0) then
                   call cerrorMsg('cannot open file', 1)
                   call cerrorMsg(input, 0)
                endif
                CALL DT_DTUINI(
     *            NEVTS,EPN,NPMASS,NPCHAR,NTMASS,NTCHAR,IDP,IEMU)
                close(TempDev)
             endif
          endif


!              This is to form  Cmsp
          call chncol(pj2, tg, a, ntp, icon)
          if( icon .ne. 0 ) then
             write(msg,*) ' cms cannot be formed' 
             call cerrorMsg(msg,0)
          endif

          cmspsave = Cmsp


!             fix model by seeing the proj. energy.
          call cfixModel(pj)
          compj = pj
          comtg = tg
          if(xbyp) then
             ek1 = pj2%fm%p(3)
          else
             ek1 = pj2%fm%p(4) -  pj2%mass
          endif
          ek = pj2%fm%p(4) -  pj2%mass

          call cbst1(1, cmspsave,  pj2, cmspj)
	  if(xbyp) then
	     ek1cms = cmspj%fm%p(3)
	  else	
	     ek1cms = cmspj%fm%p(4)-cmspj%mass
	  endif	
  	  if(pj%code .eq. kgnuc) then
            call cbst1(1, cmspsave,  pj, cmspj)
          endif	

          compjcms = cmspj
          call cbst1(2, cmspsave,  tg, cmstg)
          comtgcms = cmstg


          if(SeedFile .ne. ' ') then
             call copenfw(SeedFileDev, SeedFile, icon)
             if(icon .ne. 0) then
                call cerrorMsg(' SeedFile open error' ,0)
             endif
          endif

          write(0,*) ' ActiveModel =',ActiveMdl

          k = 1
          do while(k .le.  nevent)
             if(SeedFile .ne. ' ') then
                call rnd1s(SeedSave)
                EventNo = k
                call cwriteSeed
             endif
             if( ActiveMdl .eq. 'qgsjet2') then
!               qgsjet xs must be called
                call cxsecQGS(pj, ia, xs)
             elseif(ActiveMdl .eq. 'incdpm3') then
                call cccode2hcode(pj,  kinc)              
                xs = crossint(kinc, ek)
             endif

             if(pj%code .ne. kgnuc) then
                call chAcol(pj, ia, iz, a, ntp0)
             else
                call cheavyInt(pj, ia, iz, a, ntp0)
             endif
!/////////////
!             tgtA = ia
!             tgtZ = iz
!             call cxsecGheisha(pj, tgtA, tgtZ, xs)
!            write(0,*) xs
!
!
!            call ctotx(pj, tgtA,  xs)
!
!            write(0,*) xs
!          tgtA= 14.83   
!          call cprotonAXsec(tgtA, pj.fm.p(4)-pj.mass, xs)
!          write(0,*)sngl(sqrt(pj.fm.p(4)**2-pj.mass**2)), sngl(xs)
!////////////


!                decay sigma etc if requested
             call mydecay(a, ntp0,  ntp)

             if(Trace .ne. 0) then
                tracefile = ' '
                write(tracefile, *)
     *             TraceDir(1:klena(TraceDir))//'/trace', k
                call kseblk(tracefile, ' ', leng)
                call copenfw(TraceDev,
     *          tracefile(1:klena(tracefile)), icon)
!                   draw projectile 
                call cvisualizeTrack(compj, 1, ext)
                call cvisualizeTrack(comtg, 1, ext)
             endif

             do j = 1, ntp
                if(Trace .eq. 0 ) then
                   ptx = sngl(a(j)%fm%p(1))
                   pty = sngl(a(j)%fm%p(2))
                   pz =  a(j)%fm%p(3)
                   if(xbyp) then
                      x =  a(j)%fm%p(3) / ek1 
                   else
                      x =( a(j)%fm%p(4)- a(j)%mass)/ek1
                   endif
                   pt =sqrt(a(j)%fm%p(1)**2 + a(j)%fm%p(2)**2)
                   teta = atan2( pt, sngl(a(j)%fm%p(3)) ) *
     *               180./3.1415
                   ke = a(j)%fm%p(4)-a(j)%mass

!                    transform to CMS.
                   labp = a(j)
                   call cbst1(j, cmspsave,  a(j), a(j))
                   call cyeta(a(j), y, eta)

                   pzcms =  a(j)%fm%p(3)
                   if(xbyp) then
                      xcms =  pzcms / ek1cms 
                   else
                      xcms =( a(j)%fm%p(4)- a(j)%mass)/ek1cms
                      if(pzcms .lt. 0.) then
                         xcms= -xcms
                      endif
                   endif
                   ptcms =sqrt(a(j)%fm%p(1)**2 + a(j)%fm%p(2)**2)
                   tetacms = atan2( ptcms, sngl(a(j)%fm%p(3)) ) *
     *               180./3.1415
                   kecms = a(j)%fm%p(4)-a(j)%mass

                   write(*, '(7g14.4, 3i3, i6, 2g16.6, 2g14.4, i8)' )
     *               sngl(x), sngl(xcms), pz, pzcms, sngl(y), 
     *               sngl(eta),  pt, 
     *               a(j)%code, a(j)%subcode, a(j)%charge,
     *               ntp, teta, tetacms,
     *               ke, kecms, k
                else
                   call cvisualizeTrack(a(j), 0, ext*0.35)
                endif

             enddo
             if(Trace .ne. 0 ) then
                close(TraceDev)
             endif
             k = k + 1
          enddo   
          if(Trace .ne. 0) then
             call cerrorMsg(
     *      'Use gnuplot with "set para"; "splot filename w l"',
     *       1)
             call cerrorMsg(
     *      ' You can also use "slide" command in Util', 1) 
          endif
       end


      subroutine cvisualizeTrack(pp, inout, ext)
      implicit none
#include "Zptcl.h"
#include "Ztrackp.h"
      type(ptcl):: pp   ! input. a particle
      integer inout      ! input. 1.specifis  output order.
      real*4  ext       ! input. track length to be drawn 

      real*8 pabs
      real*4 x2, y2, z2
!               
      call cpxyzp(pp%fm, pabs) 
      if(pabs .eq. 0.) pabs = 1.


      x2 = pp%fm%p(1)/pabs*ext
      y2 = pp%fm%p(2)/pabs*ext
      z2 = pp%fm%p(3)/pabs*ext

      if(inout .eq. 0) then
         write(TraceDev, *)
     *     0., 0., 0., pp%code, sngl(pp%fm%p(4)), pp%charge
         write(TraceDev, *) 
     *     x2, y2, z2, pp%code, sngl(pp%fm%p(4)), pp%charge
      else
         write(TraceDev, *) 
     *     -x2, -y2, -z2, pp%code, sngl(pp%fm%p(4)), pp%charge
         write(TraceDev, *)
     *     0., 0., 0., pp%code, sngl(pp%fm%p(4)), pp%charge
      endif
      write(TraceDev, *)
      write(TraceDev, *)
      end

      subroutine mydecay(a, nin, n)
      implicit none
#include "Zcode.h"
#include "Zprivate.h"
#include "Zptcl.h"

!
!       particles such as sigma etc in 'a' are decayed if requested
!      
      type(ptcl):: a(nmax)  ! input/output. 

      integer nin      ! input ptcls in a.
      integer n        ! output. ptcls in a. n >= nin.

      type(ptcl):: b(nmax)   ! working array

      integer i, j, code, m
!

      m = 0
      n = nin
      do while (n .ne. m)
         m = n
         n = 0 
         do i = 1, m
            code = a(i)%code
            if( Jdecay(code) .ne. 0 ) then
               if( code .eq. kpion ) then
                  if( a(i)%charge .eq. 0 ) then 
                     call cpi0Decay( a(i),  b(n+1), j)
                  else
                     b(n+1) = a(i)
                     j = 1
                  endif
               else
                  if( code .eq. kdmes ) then
                     call cdDecay( a(i), b(n+1), j)
                  elseif( code .eq. keta ) then
                     call cetaDecay( a(i), b(n+1), j)
                  elseif( code .eq. kgzai ) then
                     call cgzaiDecay( a(i), b(n+1), j ) 
                  elseif( code .eq. klambda ) then
                     call clambdaDcy( a(i), b(n+1), j ) 
                  elseif( code .eq. klambdac ) then
                     call clambdacDcy( a(i), b(n+1), j ) 
                  elseif( code .eq. ksigma ) then
                     call csigmaDecay( a(i), b(n+1), j ) 
                  elseif( code .eq. kbomega ) then
                     call cbomegaDcy( a(i), b(n+1), j ) 
                  else
                     write(*,*)  ' code =', code
                     call cerrorMsg('mydecay error', 0)
                  endif
               endif
               n = n + j
            else
               b( n+1 )  = a(i)
               n = n + 1
            endif 
            if( n .gt. nmax ) then
               call cerrorMsg('partilce array overflow', 0)
            endif
         enddo
         if(n .ne. m) then
            do i = 1, n
               a(i) = b(i)
            enddo
         endif
      enddo
      end

