#include "ZepMaxdef.h"
#include "Zelement.h"
#include "ZbpTbl.h"
#include "Zstern.h"         
#include "ZbpSample.h"
#include "ZbpPE.h"
#include "Zurban.h"
#include "Zmubpn.h"
#include "Zxcom.h"


       integer maxElements
       parameter (maxElements=MAX_ELEMENTS)
       type epmedia 
       sequence

       integer noOfElem          ! actual number of elements
       type(element):: elem(maxElements)
       real*8  No(maxElements)  ! number of each element; copied to
               !  OrigNo and then normalized. so that sum be 1.0
       real*8  OrigNo(maxElements)  ! number of each element            
!     ***   once tried to change this to be contained in elem and
       
!     use like media%elem(i)%No (seems more logical)  than media%No(i)
!     and consistent with media%elem(i)%A etc.       
!           but keep the old style.  2018/Oct.5 ***
!        in the molecule; if a mixture of
               ! molecules/atoms,(volume fraction of each molecule/atom) times 
               ! number of each element  in the molecule.
               ! E.g:  H2O --> 2, 1 .  Air, say,  79% N2, 20% O2 and 1% Ar
               ! by volume, 0.79x2, 0.20x2, 0.01
               ! This is later normalized so that the sum be 1.0
       real*8  w(maxElements)   !  No(i)A(i)/sum(No(i)A(i))  same note as No
               ! next is added from v8.77
       real*8  npercm3(maxElements) ! # of i-th element /cm3=No/Ai*rho*wi same
                      ! note as No
       real*8  nsigma(maxElements) ! No(i)s_i
       real*8  sumns            !  sum of above

       real(8):: sumNo   ! sum of OrigNo(:) 				 
!!       real*8  sigma(nxsec)     !  sum No(j)*elem(j).sigma(i), i=1,nxsec
!c            real*8  cumsigma(maxElements, nxsec) !
!c                 !      cumulative  sum(j) of 
!c                 !      No(j).elem(j).sigma(i)/sigma(i) for j=1, nxsec
!c                 !
       integer colElem          ! element # at which interaction took place
       integer colA             !  A of such one:. int(elem(colElm)%A+0.5)
       integer colZ             !  Z of such one:  int(elem(colElm)%Z)
       real(8):: colXs            ! x-section for that target(mb)
       real(8):: xs   ! this is xs for the media.
       
            real*8  ndensity         !  effective number density /cm^3
            real*8  wp    ! plasma frequency  x hbar (GeV)
            real*8  n     ! refractive index
            real*8  nd    ! number of ingredients / cm^3 
            real*8  A     ! sum No x  Ai
            real*8  Z     ! sum No x  Zi
            real*8  Z2    ! sum No x  Zi**2
            real*8  ZZ1   ! sum No x  Zi(Zi+1)
               ! for electron;  this * t /(gamma beta2)^2 = Xc^2 (t g/cm2)
            real*8  MoliereForXc2
            real*8  MoliereExpb !  exp(b) = t x this  (t in g/cm2)a
                  ! for z=1 and beta =1.          
                  ! this  = 6702 sum/A
                  !    sum = Sum No x Zi^(1/3)(Zi+1)/(1+3.327(Ziz/137)^2)
            real*8  Z1_3rd !  <Z^1/3>  not <Z>^(1/3)
            real*8  Z2_3rd !  <Z^2/3>  not <Z>^(2/3) 
            real*8  mbtoPgrm ! 10^-27 x N0/A. If multiplied to sigma in mb,
                           !  we obtain probability / (g/cm2).
            real(8):: mbtoPkgrm ! mbtoPgrm/10d0
            real*8  mbtoPcm  ! rho x mbtoPgrm. If multiplied to sigma in mb, 
                           ! we obtaind probability / cm
            real*8  mbtoPX0  ! mbtoPgrm x X0g.  If multiplied to sigma in mb, 
                             !  we obtain probability /radation length.
!                            next ones are used when we approximate
!               a compound /molecule as an atom
            real*8  mbtoPgrm2
            real*8  mbtoPcm2
            real*8  mbtoPX02 

            real*8  Z2byAeff !  sum wi x  Zi**2/Ai
            real*8  Z5byAeff !  sum wi x  Zi**5/Ai
            real*8  Aeff  ! sum wi x Ai
            real*8  Z2eff ! Z2byAeff x Aeff
            real*8  Zeff  ! sqrt(Z2eff)
            real*8  Zeff3 ! Zeff**(1/3)
            real*8  LogZ  ! log(Zeff) 
            real*8  A2eff ! sum wi x Ai^2
            real*8  ZbyAeff ! sum wi x Zi/Ai
            real*8  I     !  average ionization potential energy in GeV.
            real*8  rho   ! density in g/cm^3
            real*8  X0    ! radiation length.  in cm
            real(8):: X0m ! radiation lenght in m.
            real*8  X0g   ! radiation length.  in g/cm^2
            real(8)::X0kg ! reaiation length in  kg/m^2 =X0g*10
            real*8  gtocm ! g/cm^2 to cm.
	    real(8)::kgtom !  gtocm*1.0d-3: kg/m2 to m.			  
            real*8  dEdxatp3m ! dE/dx at p=3me for electron. ~ Ecrit
            real(8):: Ecrit     !  GeV  for electron
            real(8):: Ecritmu    !  GeV  for muon
            real*8  rhoc  ! comp.rhoc is copied whenever new comp. comes
 		   ! note;this is real*8 while comp.rhoc is real*4
            integer gasF  ! flag for gas. If  1, media is gas,  0 -->solid
            character*8  name  ! name of media
            integer  format    ! format of the  basic table. (1 or 2)

            real*8  s1         ! Migdal's s1
            real*8  logs1      ! log(s1) 
            real*8  basearea   ! pi x Re**2 * N* Z /A *X0g  = 0.15 Z/A*X0g
            real*8  cScrC1     ! const which appears in the complete screening
                               ! crossection
            real*8  cScrC2     ! the other such one
            real*8  cScrMain   ! (4/3C1 + C2)
            real*8  BirksC1    !  quenching correction coef.
            real*8  BirksC2    !  //
            real*8  BirksCC    !  //
            character(1)::Birks ! flag to identify what quenching correction
                                ! should be applied using BirksC1, etc.
              !  In the media  data,  data format  for c1, c2, cc is one of 
              !  a)  c1 c2 cc   (none zero or all zero)
              !  b)  c1 c2 cc  B
              !  c)  c1 c2  T
              !  d)  c1 c2 cc  L 
              !  a) is old format. all  0 means no quencthing.
              !     if not so, or b) format case,  only c1 will be used
              !     as the original formula by Birks
              !         dE/dx|q =  dE/dx/(1+c1*dE/dx)
              !     e.g. c1 (13 g/cm2/GeV) for organic scinit 
              !          c2 (9.6 (g/GeV/cm2)^2 .. is intended to use
              !       (1 + c1 X + c2*X**2) as the denominator.
              !       where X is dE/dx in GeV/(g/cm2). 
              !      cc is the correction factor for heavy particles
              !       (e.g 0.5714).  c2 and cc are not used in Epics
              !  c)  Tarle's formula  dE/dx|q = (1-c2)X/( 1+c1(1-c2)X) + c2X
              !              c1 in (g/cm2)/GeV. c2 simple number
              !         if fitting is done in another unit, simply change c1
              !         (say, if X is MeV/(g/cm2), C1*1000 -->c1)
              !  d)  Log formula
              !      dE/dx|q = X * (cc*X)**(-c2*log(c1*X))
              !      c1, cc  in g/cm2/GeV c2 simple number.
              !      if fitting is done in another unit, change 
              !      c1 and cc. (e.g if in MeV/(g/cm2). C1*1000->c1
              !      CC*1000->cc
            integer srim   !  index for srim data in module srimdata
       type(bpTbl)::   tbl
       type(sternh)::  sh
       type(SmpCnst)::  cnst
       type(photoE)::  pe
       type(urban)::  urb
       type(mubpn)::  mu 
       type(epxcom)::  xcom
       end type epmedia
