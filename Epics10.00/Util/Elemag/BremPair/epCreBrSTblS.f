!   #define DEBUGS
      subroutine epCreBrSTblS(mediain, cnst)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"
!     *                                                             
!     *       create sampling table for Bremstrahlung 
!     *       at low energies using Seltzer's data

!
       type(epmedia):: mediain  !  media
       type(SmpCnst)::  cnst  ! must be media.cnst
!
!        cp mediain into common area
      media = mediain
      call epwtBrCnstS(cnst)
      call epCreBrSeltz1(cnst)
      call epCreBrSeltz2(cnst)
      end  subroutine epCreBrSTblS



      subroutine epCreBrSeltz1(cnst)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"
       type(SmpCnst)::  cnst  ! must be media.cnst
!            lower energy region 
!
      integer ws
      parameter (ws = 10000)
      real*8 work(ws), work2(ws)

!            total x-sec.
      if( cnst%BrTXTS  .gt. mxBrTXS) then
         write(0,*) 'cnst%BrTXTS=',cnst%BrTXTS,
     *  "> mxBrTXS=",  mxBrTXS
         call cerrorMsg(
     *   'too large total X-sec. table for lower Seltzer brems', 0)
      endif
      call epCreBrTXTS(cnst, work, work2, cnst%BrTXTS)

!           sampling table in A region
      if( cnst%BrUszSA * cnst%BrES .gt. mxBrTblSA .or.
     *    cnst%BrUszSA * cnst%BrES .gt. ws ) then
         write(0,*)
     *  "cnst%BrUszSA * cnst%BrES =", cnst%BrUszSA,"*",
     *  cnst%BrES, "=", cnst%BrUszSA*cnst%BrES,
     *  ">  mxBrTblSA=", mxBrTblSA
         write(0,*) ' ws =',ws
         call cerrorMsg('brem table for SA is too large', 0)
      endif
      call  epCreBrSA(cnst, work, cnst%BrUszSA, cnst%BrES)
!           sampling table in B region
      if( cnst%BrUszSB * cnst%BrES  .gt. mxBrTblSB) then
         write(0,*)
     *    "cnst%BrUszSB * cnst%BrES=",
     *     cnst%BrUszSB, "*",  cnst%BrES, " =",
     *    cnst%BrUszSB * cnst%BrES, "> mxBrTblSB=",  mxBrTblSB
         call cerrorMsg('brem table for SB is too large', 0)
      endif
      call  epCreBrSB(cnst, work, cnst%BrUszSB, cnst%BrES)
      end  subroutine epCreBrSeltz1

      subroutine epCreBrSeltz2(cnst)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"
       type(SmpCnst)::  cnst  ! must be media.cnst
!            higher energy region 
!
      integer ws
      parameter (ws = 10000)
      real*8 work(ws), work2(ws)

!            total x-sec.
      if( cnst%BrTXTS2  .gt. mxBrTXS2) then
         write(0,*) 'cnst%BrTXTS2=',cnst%BrTXTS2,
     *  "> mxBrTXS2=",  mxBrTXS2
         call cerrorMsg(
     *   'too large total X-sec. table for higher Seltzer brems', 0)
      endif
      call epCreBrTXTS2(cnst, work, work2, cnst%BrTXTS2)

!           sampling table in A region
      if( cnst%BrUszSA2 * cnst%BrES2 .gt. mxBrTblSA2) then
         write(0,*)
     *  "cnst%BrUszSA2 * cnst%BrES2=", cnst%BrUszSA2,"*",
     *  cnst%BrES2, "=", cnst%BrUszSA2*cnst%BrES2,
     *  ">  mxBrTblSA2=", mxBrTblSA2
         call cerrorMsg('brem table for SA2 is too large', 0)
      endif
      call  epCreBrSA2(cnst, work, cnst%BrUszSA2, cnst%BrES2)
!           sampling table in B region
      if( cnst%BrUszSB2 * cnst%BrES2  .gt. mxBrTblSB2) then
         write(0,*)
     *    "cnst%BrUszSB2* cnst%BrES2=",
     *     cnst%BrUszSB2, "*",  cnst%BrES2, " =",
     *    cnst%BrUszSB2 * cnst%BrES2, "> mxBrTblSB2=", mxBrTblSB2
         call cerrorMsg('brem table for SB2 is too large', 0)
      endif
      call  epCreBrSB2(cnst, work, cnst%BrUszSB2, cnst%BrES2)
      end  subroutine epCreBrSeltz2

!     ****************************************
      subroutine epCreBrTXTS(cnst,  erg, tbl, size)
!     ****************************************
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer size
      real*8  tbl(size), erg(size)

      real*8  Ek, Ee,  vmax, vmin,  tcb, de
      character*160 msg
      integer i

!   /////////////////////////
!    continuation at 100 MeV Seltzer and partial screening
!    formula. generally heavy materials are good
!    low Z material has glitch of 10 -20  %
!    At present no normalization is tried.
!    To get the value below, uncomment the
!    part of the  program below.
!==================== at 1GeV
!Pb
! Brem XS at E=  0.999990000000000        614937.362267548     
! Brem XS at E=   1.00001000000000        624285.473828595    
!
!Al
!  Brem XS at E=  0.999990000000000        21685.4116313612     
!  Brem XS at E=   1.00001000000000        21679.3370821324   
!BGO
!  Brem XS at E=  0.999990000000000        157300.388762514     
!  Brem XS at E=   1.00001000000000        159377.105885142    
!
!Si
!  Brem XS at E=  0.999990000000000        24865.9911773002     
!  Brem XS at E=   1.00001000000000        24834.9687550972    
!Air
!  Brem XS at E=  0.999990000000000        7770.07544354533     
!  Brem XS at E=   1.00001000000000        7665.26824930834
!
! plastic Scinti
!  Brem XS at E=  0.999990000000000        2924.89850986149     
!  Brem XS at E=   1.00001000000000        2860.61100050295 
!
! Fe
!  Brem XS at E=  0.999990000000000        78256.5221115459     
!  Brem XS at E=   1.00001000000000        78104.3555432043 
! W
!  Brem XS at E=  0.999990000000000        509197.902670278     
!  Brem XS at E=   1.00001000000000        515653.327859797 
! Cu
!  Brem XS at E=  0.999990000000000        95253.1157581748     
!  Brem XS at E=   1.00001000000000        95604.5329845806
! H2 gas
!  Brem XS at E=  0.999990000000000        335.977047803700     
!  Brem XS at E=   1.00001000000000        306.555938740378
!
! ===================== at 100 MeV
! BGO
! Brem XS at E=  9.999900000000000E-002  114622.885437669
! Brem XS at E=  .100001000000000  115022.240288808
!
!Pb
! Brem XS at E=  9.999900000000000E-002  452171.630320046
! Brem XS at E=  .100001000000000  455365.448006682
!
!H2 (gas)
! Brem XS at E=  9.999900000000000E-002  268.026568486841
! Brem XS at E=  .100001000000000  217.605740594913
!
!Si
! Brem XS at E=  9.999900000000000E-002  18296.4484918632
! Brem XS at E=  .100001000000000  17838.2578766329
!
!Fe
! Brem XS at E=  9.999900000000000E-002  57120.4530869126
! Brem XS at E=  .100001000000000  56147.6825177278
!
!W
! Brem XS at E=  9.999900000000000E-002  378705.118797284
! Brem XS at E=  .100001000000000  380306.005877633
!
!Al
! Brem XS at E=  9.999900000000000E-002  15981.7824199197
! Brem XS at E=  .100001000000000  15567.2559021630
!
!pla Scinti
! Brem XS at E=  9.999900000000000E-002  2197.66802947518
! Brem XS at E=  .100001000000000  2043.99285184687
!
!Cu
! Brem XS at E=  9.999900000000000E-002  69520.4358663770
! Brem XS at E=  .100001000000000  68777.3779982839
!  ===========================
!     continuation at 10 MeV case
! BGO
!  Brem XS at E=   9.99990000000000D-003   67272.8642710921     
!  Brem XS at E=   1.00001000000000D-002   65514.8334236348   
! H2(gas)
!  Brem XS at E=   9.99990000000000D-003   124.067255033441     
!  Brem XS at E=   1.00001000000000D-002   115.381548133051  
! Si
!  Brem XS at E=   9.99990000000000D-003   10152.6381914789     
!  Brem XS at E=   1.00001000000000D-002   10003.6693393773   
! Al  
!  Brem XS at E=   9.99990000000000D-003   8834.05295739954     
!  Brem XS at E=   1.00001000000000D-002   8717.92352096671   
! Cu
!  Brem XS at E=   9.99990000000000D-003   39683.8031823851     
!  Brem XS at E=   1.00001000000000D-002   39003.4735041571    
! plasticSCIN
!  Brem XS at E=   9.99990000000000D-003   1151.25081712306     
!  Brem XS at E=   1.00001000000000D-002   1121.94237263980 
! W 
!  Brem XS at E=   9.99990000000000D-003   222687.620819852     
!  Brem XS at E=   1.00001000000000D-002   217007.840876563   
! Pb  
!  Brem XS at E=   9.99990000000000D-003   267092.652128783   
!  Brem XS at E=   1.00001000000000D-002   259782.820635386   
!

      write(msg,*) ' Matter=', media%name,
     *  ': Computing the total Xsec of Seltzer Brems from Ee= ',
     *    cnst%BrEeminS, ' GeV'
      call cerrorMsg(msg, 1)
      Ek = cnst%BrEeminS - masele  !  Ek

      de = 10.d0**cnst%BrdETXS
      do   i=1, size
         Ee = Ek + masele
         call epChkEndValSel(Ee, i, size)
         Eeme = Ee/masele
!         vmax=1.d0
!         vmin=cnst.BrEgminS/(E-masele)
         vmax = 1.d0 - masele/Ee   !  Eg/Ee= (Ee-masele)/Ee = 1-masele/Ee
         vmin= media%cnst%BrEgminS/Ee  !  BrEgminS is not ratio 
                                    ! in v9.135
         if(vmin .ge.  vmax) then
            call cerrorMsg('Eg/E <= vmin in lower Seltzer tab', 0)
         endif
         call epBrgeneTX(vmin,  vmax,  tcb)
         tbl(i) = tcb* media%mbtoPX0  ! prob. per radiation length.
         erg(i) = Ek
         Ek = Ek * de
      enddo
      write(msg, *) 'Table has been made up to  Ek=', Ek/de,' GeV'
      call cerrorMsg(msg, 1)
      write(msg, *) 'Seltzer Brems total X-sec. table upto Ek=',
     *     Ek/de,' GeV'
      call epwt1dTbl(msg, erg,  tbl, size, media%name)

      end
!     ****************************************
      subroutine epCreBrTXTS2(cnst,  erg, tbl, size)
!     ****************************************
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer size
      real*8  tbl(size), erg(size)

      real*8  Ek, Ee,  vmax, vmin,  tcb, de
      character*160 msg
      integer i

      write(msg,*) ' Matter=', media%name,
     *  ': Computing the total Xsec of Seltzer Brems from Ee= ',
     *    cnst%BrEeminS2, ' GeV'
      call cerrorMsg(msg, 1)
      Ek = cnst%BrEeminS2 - masele  !  Ek

      de = 10.d0**cnst%BrdETXS2
      do   i=1, size
         Ee = Ek + masele
         call epChkEndValSel2(Ee, i, size)
         Eeme = Ee/masele
!         vmax=1.d0
!         vmin=cnst.BrEgminS/(E-masele)
         vmax = 1.d0 - masele/Ee   !  Eg/Ee= (Ee-masele)/Ee = 1-masele/Ee
         vmin= media%cnst%BrEgminS2  !  BrEgminS2 is ratio e.g 10^-5

         if(vmin .ge.  vmax) then
            call cerrorMsg('Eg/E <= vmin in higher  Seltzer tab', 0)
         endif
         call epBrgeneTX(vmin,  vmax,  tcb)
         tbl(i) = tcb* media%mbtoPX0  ! prob. per radiation length.
         erg(i) = Ek
         Ek = Ek * de
      enddo
      write(msg, *) 'Table has been made up to  Ek=', Ek/de,' GeV'
      call cerrorMsg(msg, 1)
      write(msg, *) 'Seltzer Brems total X-sec. table upto Ek=',
     *     Ek/de,' GeV'
      call epwt1dTbl(msg, erg,  tbl, size, media%name)

      end
!     **********************************************
      subroutine epCreBrSA(cnst, bla, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee 
      real*8  bla(sizeu,  sizee)

      real*8 de1, Ek, Ee, vmin, vl, vr, u, tcb, vx, tcbx,
     *  eps, v

      external epBrgeneSolv
      real*8  epBrgeneSolv

      real*8  epBrgene

      character*160 msg
      integer i, iu, j


      common/upsic/upsi, vmax
      real*8 vmax, upsi

      data eps/1.d-8/


      de1 = 10.d0**cnst%BrdES

      Ee = cnst%BrEeminS    
      Ek = Ee - masele
      write(msg, *)
     *  ' Creating lower Seltzer Brem sampling table: Ek>=',
     *    Ek
      call cerrorMsg(msg, 1)
      call cerrorMsg('for small Eg region',1)

      do  i = 1, sizee
         Ee = Ek + masele
!           chech env value
         call epChkEndValSel(Ee, i, sizee)

         Eeme = Ee/masele
!         vmin=cnst.BrEgminS/(E-masele)
!         vmax=1.d0
!         vmin=cnst.BrEgminS/E
         vmin= media%cnst%BrEgminS/Ee ! now EgminS is not ratio
         vmax=1.d0 - masele/Ee
         call epBrgeneTX(vmin, vmax, tcb)
         vl=vmin
         vr=vmax
         
         u=cnst%BrUminSA
         
         do  iu=1, sizeu -1
            
            upsi=u*tcb
!             v here obtained is Eg/Ee
            call kbchop(epBrgeneSolv, vl, vr, eps, v, j)
            if(j .le. 0) then
               
               write(0, *) ' cond, Ek, u=', j, Ek, u
               do vx = vmin, vmax, (vmax-vmin)/100.d0
!                  call eptotcbS(vx, vmax, tcbx)
                  call epBrgeneTX(vx, vmax, tcbx)
                  write(0, *) vx, tcbx
               enddo
               write(0,*) 'u upsi, vl, vr', u,
     *              upsi, vl, vr,
     *              ' vmin=', vmin, ' tcb=', tcb
            endif
            bla(iu,i)= log(v/vmin)/(1.-u)
#if defined (DEBUGS)            
!!!!!!
            write(0,'(a, i3,1p, 4g14.6)')
     *           'ss1 ', iu, Ek, u, v, bla(iu, i)
!!!!!!!            
#endif
            u=u + cnst%BrdUSA
            vr = v
         enddo
         bla(sizeu, i)=tcb/vmin/epBrgene(media, force, Ee, vmin)
#if defined (DEBUGS)                     
!!!!!!!!
         write(0,'(a, i3,1p, 4g14.6)')
     *           'ss1 ', sizeu, Ek, u, vmin, bla(sizeu, i)
!!!!!!!
#endif         
         Ek =Ek* de1
      enddo
      call epwt2dTbl(
     *  'Lower Selzter Brems sampling table; log(v/vmin)/(1-u)',
     *   bla, sizeu, sizee)
      call cerrorMsg('the table has been created', 1)
      end
!     **********************************************
      subroutine epCreBrSA2(cnst, bla, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee 
      real*8  bla(sizeu,  sizee)

      real*8 de1, Ek, Ee, vmin, vl, vr, u, tcb, vx, tcbx,
     *  eps, v

      external epBrgeneSolv
      real*8  epBrgeneSolv

      real*8  epBrgene

      real(8),parameter:: dv=10.0d0**0.05d0
      character*160 msg
      integer i, iu, j


      common/upsic/upsi, vmax
      real*8 vmax, upsi

      data eps/1.d-8/
      de1 = 10.d0**cnst%BrdES2

      Ee = cnst%BrEeminS2    
      Ek = Ee - masele
      write(msg, *)
     *  ' Creating higher Seltzer Brem sampling table: Ek>=',
     *    Ek
      call cerrorMsg(msg, 1)
      call cerrorMsg('for small Eg region',1)

      do  i = 1, sizee
         Ee = Ek + masele
!           chech env value
         call epChkEndValSel2(Ee, i, sizee)

         Eeme = Ee/masele
!         vmin=cnst.BrEgminS/(E-masele)
!         vmax=1.d0
!         vmin=cnst.BrEgminS/E
         vmin= media%cnst%BrEgminS2 !  EgminS2 is  ratio
         vmax=1.d0 - masele/Ee
         
         call epBrgeneTX(vmin, vmax, tcb)
         vl=vmin
         vr=vmax
         
         u=cnst%BrUminSA2
         
         do  iu=1, sizeu -1
            
            upsi=u*tcb
!             v here obtained is Eg/Ee
            call kbchop(epBrgeneSolv, vl, vr, eps, v, j)
            if(j .le. 0) then
               
               write(0, *) ' cond, Ek, u=', j, Ek, u
               write(0,*) 'u upsi, vl, vr', u,
     *              upsi, vl, vr,
     *              ' vmin=', vmin, ' tcb=', tcb
                              
               vx = vmin
!     do vx = vmin, vmax, (vmax-vmin)/100.d0
               do while ( vx < vmax)
!                  call eptotcbS(vx, vmax, tcbx)
                  call epBrgeneTX(vx, vmax, tcbx)
                  write(0, *) vx, tcbx, 1.d0 - tcbx/upsi
                  vx = vx * dv
               enddo
            endif
            bla(iu,i)= log(v/vmin)/(1.-u)
#if defined (DEBUGS)                                 
!!!!! 
            write(0,'(a, i3, 1p, 4g13.4)') 'ss2 ', iu, Ek, u, v,
     *      bla(iu,i)
!!!!!!!!!!
#endif            
            u=u + cnst%BrdUSA2
            vr = min(v*dv*dv*dv, vmax)
         enddo
!         bla(sizeu, i)=tcb/vmin/epBrSfs(media, Eeme, vmin) ! <= v9.136
         bla(sizeu, i)=tcb/vmin/epBrgene(media, force, Ee, vmin)
#if defined (DEBUGS)                                          
!!!!!!!!!         
         write(0,'(a, i3, 1p, 4g13.4)') 'ss2 ', sizeu, Ek, u, vmin,
     *        bla(sizeu,i)
!!!!!!!!
#endif         
                           ! generic  function must be used otherwise, when LPM
                           ! comes in this region, results in wrong.
                           ! (eBrSfs dose not include LPM)

         Ek =Ek* de1
      enddo

      call epwt2dTbl(
     *  'Higher Selzter Brems sampling table; log(v/vmin)/(1-u)',
     *   bla, sizeu, sizee)
      call cerrorMsg('the table has been created', 1)
      end

!     ***********************************************
      subroutine epCreBrSB(cnst, blb, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee 
      real*8  blb(sizeu, sizee)

      real*8 de1, Ek, Ee, vmin, vl, vr, u, tcb, 
     *   eps, us, v

      external epBrgeneSolv
      real*8  epBrgeneSolv
      real*8  epBrgene

      integer  iu, j, ie


      common/upsic/upsi, vmax
      real*8 vmax, upsi


      call cerrorMsg(
     * 'Creating brems sampling table by Seltzer data for',1)
      call cerrorMsg(' large Eg ',1)


      Ee = cnst%BrEeminS 
      Ek = Ee - masele
      de1 = 10.d0**cnst%BrdES
      do ie = 1, sizee
         Ee = Ek + masele

         call epChkEndValSel(Ee, ie, sizee)

         Eeme = Ee/masele
!          vmin=cnst.BrEgminS/(E-masele)
!          vmax=1.d0
!          vmin=cnst.BrEgminS/E
         vmin=cnst%BrEgminS/Ee  ! EgminS is not ratio in v9.135
         vmax=1.d0- masele/Ee
         call epBrgeneTX(vmin, vmax, tcb)
         vl=vmin
         vr=vmax
         us=cnst%BrUminSB + cnst%BrdUSB
         do   iu = 2,  sizeu
            u=us**4
            upsi=u*tcb
            if(iu .le .6) then
               eps = 1.d-8
            else
               eps = 1.d-8
            endif
            
            call kbchop(epBrgeneSolv, vl, vr, eps, v, j)
            if(j .le. 0) then
               write(0,'('' Ek,u='',2g12.3)') Ek, u
            endif
            blb(iu,ie) = log(vmax/v)/u
#if defined (DEBUGS)                                                      
!!!!! 
            write(0,'(a, i3, 1p, 5g13.4)') 'ss3 ', iu, Ek, us, v,
     *             blb(iu,ie), u
!!!!!!!!!!
#endif            
            us=us + cnst%BrdUSB
            vr = v
         enddo
         blb(1, ie) = tcb/vmax/epBrgene(media, force, Ee, vmax)
#if defined (DEBUGS)                                                      
!!!!! 
         write(0,'(a, i3, 1p, 5g13.4)') 'ss3 ', 1, Ek, cnst%BrUminSB,
     *     vmax,  blb(1, ie), cnst%BrUminSB**4
!!!!!!!!!!
#endif                     
         Ek= Ek*de1
       enddo
       
       call epwt2dTbl(
     *  'Lower Seltzer Brems sampling tbl in B; log(v/vmin)',
     *  blb, sizeu, sizee)
       call cerrorMsg('the table has been created', 1)
       end
!     ***********************************************
      subroutine epCreBrSB2(cnst, blb, sizeu,  sizee)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

       type(SmpCnst)::  cnst  ! must be media.cnst
      integer sizeu, sizee 
      real*8  blb(sizeu, sizee)

      real*8 de1, Ek, Ee, vmin, vl, vr, u, tcb, 
     *   eps, us, v

      external epBrgeneSolv
      real*8  epBrgeneSolv
      real*8  epBrgene

      integer  iu, j, ie


      common/upsic/upsi, vmax
      real*8 vmax, upsi


      call cerrorMsg(
     * 'Creating brems sampling table by higher Seltzer data for',1)
      call cerrorMsg(' large Eg ',1)


      Ee = cnst%BrEeminS2 
      Ek = Ee - masele
      de1 = 10.d0**cnst%BrdES2
      do ie = 1, sizee
         Ee = Ek + masele

         call epChkEndValSel2(Ee, ie, sizee)

         Eeme = Ee/masele
         vmin=cnst%BrEgminS2  ! EgminS2 is ratio
         vmax=1.d0- masele/Ee
         call epBrgeneTX(vmin, vmax, tcb)
         vl=vmin
         vr=vmax
         us=cnst%BrUminSB2 + cnst%BrdUSB2
         do   iu = 2,  sizeu
            u=us**4
            upsi=u*tcb
            if(iu .le .6) then
               eps = 1.d-8
            else
               eps = 1.d-8
            endif

            call kbchop(epBrgeneSolv, vl, vr, eps, v, j)
            if(j .le. 0) then
               write(0,'('' Ek,u='',2g12.3)') Ek, u
            endif
            blb(iu,ie) = log(vmax/v)/u
#if defined (DEBUGS)
!!!!!
            write(0,'(a, i3, 1p, 5g13.4)') 'ss4 ', iu, Ek, us, v,
     *           blb(iu,ie), u
!!!!!!!!!!            
#endif
            us=us + cnst%BrdUSB2
            vr = v
         enddo
         blb(1, ie) = tcb/vmax/epBrgene(media, force, Ee, vmax)
#if defined (DEBUGS)                                                           
!!!!! 
         write(0,'(a, i3, 1p, 5g13.4)') 'ss4 ', 1, Ek, cnst%BrUminSB2,
     *         vmax,    blb(1, ie), cnst%BrUminSB2**4
!!!!!!!!!!
#endif
         Ek= Ek*de1
       enddo
       
       call epwt2dTbl(
     *  'Higher Seltzer Brems sampling tbl in B; log(v/vmin)',
     *  blb, sizeu, sizee)
       call cerrorMsg('the table has been created', 1)
       end

      subroutine  epChkEndValSel(Ee, i, size)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

      real*8 Ee  ! in/out
      integer i, size

!            table end value should not exceed the region end value
      if(i .eq. size) then 
         if(Ee .ge. media%cnst%BrEemaxS) then
            if(Ee .gt. media%cnst%BrEemaxS*1.0001d0) then
               write(0,*) '*********** too large end value *** '
               write(0,*) 'Ee=',Ee, 
     *           '  > media%cnst%BrEemaxS*1.001d0=',
     *           media%cnst%BrEemaxS*1.0001d0
               stop
            endif
            Ee = media%cnst%BrEemaxS*0.999999999d0
         endif
      endif
      end
      subroutine  epChkEndValSel2(Ee, i, size)
      implicit none
#include "Zmedia.h"
#include "Zmass.h"
#include "ZBPgene.h"

      real*8 Ee  ! in/out
      integer i, size

!            table end value should not exceed the region end value
      if(i .eq. size) then 
         if(Ee .ge. media%cnst%BrEemaxS2) then
            if(Ee .gt. media%cnst%BrEemaxS2*1.0001d0) then
               write(0,*) '*********** too large end value *** '
               write(0,*) 'Ee=',Ee, 
     *           '  > media%cnst%BrEemaxS*1.001d0=',
     *           media%cnst%BrEemaxS2*1.0001d0
               stop
            endif
            Ee = media%cnst%BrEemaxS2*0.999999999d0
         endif
      endif
      end
