#include "ZepMaxdef.h"
  ! if you change the Earth to another world, cputCerenkov.f in
  ! Tracking might be updated, if you use it for Cerenkov light
  ! generation.
module modAtmosDef
  implicit none
  !     this should be linked to Cosmos/Moddule/catmosDef.f90
  !     as ln -s ../Atmoshere/catmosDef.f90 catmosDef.f90
  !     in Cosmos/Module/
  !  This module defines basic variables to describe atmosphere
  !  and can be used for non-Earth atmosphere too.
  !  Many variables in Zearth.h and Zatmos.h in Cosmos/cosmos
  !  were  moved and updated.
  !
  !
  
  real(8),save:: Eradius = 6378136.0d0  !  radius in m
!  real(8),parameter:: Eecen2= 0.006694470d0  !  square of the ecentricity
  ! not used at present
!  real(8),save:: Eone_ecen2 = 1.d0 - Eecen2 !   1 - ecen2
!  real(8),save:: Maxheight=220.d3 !  max height with some air. m  OBSO
  real(8),save:: Minheight = -500.d0 !min height without soil. m

  
  real(8),save:: MagChgDist =20.d3   !2 Distance where mag. can be seen
  ! as const.(m) at sea level (< 1% diff.)  Moved from
  ! modEMcontrol in cepControl.f Jul.27,2019
  real(8),save:: gaccel = 9.80665d0   ! gravitational accel. m/sec2
                  ! at near the surface of the Earth
  !----------------------- from Zatmos.h >>>>>>>>>>>
  real(8),save:: Hinf = 200.d3 ! m  At very high altitude
  ! (> Maxheight),  air density is calulated using
  ! this value as the scale height.  rho= rho_m exp( -(z-Zm)/Hinf)
  !     -M means the value at Maxheight.
    !
  real(8),save:: AlmostVacH = 210d3   ! m
  real(8),save:: AlmostVacT   ! = 1.8d-3  !  kg/m2  will be set in
                       !  catmosCnst1 using AlmostvacH
  !
  integer,parameter:: maxnodes = 40
  ! maximum no. of nodes to express the atmosphere.
  !      We inverse the following relation to be
  !  applicable to non-Earth case.
  !       for length to thickness conversion and  v.v 
!  real(8),parameter:: LenStep=1000.d0   ! in m
!  integer,save:: maxl2t=1200.0d3/LenStep !  = 1200
     ! slant length 1200km ~ vertical length 100km for
  !  cos~0.083.  The length is divided by LenSetp.
  integer,parameter::maxl2t=1200  !  see the bottom how this one is used
  integer:: LenStepi
  real(8),save:: LenStep
  real(8),save:: Maxheight=220.d3 !  max height with some air. m
     
!      When using the table,
!      Non-integral point (i.e., values not in the table) is
!      obtained by using interpolation using 3 integral points.
!      (see cl2tTblA.f)
!     maxnodes: 

  integer,save:: NumStep ! table size made by sub.  cl2tTbl
  real(8),save:: Htop, Hbase
  real(8),save:: Zsave, ThickSave
  real(8),allocatable,save::ThickTbl(:), LenTbl(:), CosTbl(:), &
         HeightTbl(:)  ! size should be maxl2t

  
  type atmosph 
     sequence
     real*8  z(maxnodes)
     real*8  zi(maxnodes)
     real*8  T(maxnodes)
     real*8  P(maxnodes)
     real*8  logP(maxnodes)
     real*8  rho(maxnodes)
     real*8  logrho(maxnodes)
     real*8  logrhoi(maxnodes) 
     real*8  H(maxnodes)
     real*8  a(maxnodes)
     real*8  b(maxnodes)
     real*8  cumd(maxnodes)
     real*8  logcumd(maxnodes)
     real*8  logcumdi(maxnodes)
     real*8  d0(maxnodes)
     real*8  coefh2r(maxnodes,3)
     real*8  coefh2T(maxnodes,3)
     real*8  coefh2d(maxnodes,3)
     !        real*8  coefh2P(maxnodes,3)      
     real*8  coefh2H(maxnodes,3)      
     !        real*8  coefr2h(maxnodes,3)      
     real*8  coefd2r(maxnodes,3)  
     real*8  coefd2h(maxnodes,3)      
        !        real*8  coefP2h(maxnodes,3)
     integer:: node2mediaNo(maxnodes) !  for a given node number, n,
        !  j= node2mediaNo(n) gives index for media (for usual atmphere
        !  1)  and media(j) gives the media name (normally "air")
        !  max of j can be nodes<=maxnodes.
        !  If, for a give n, "air" is specified, it means meidia
        !  between node n and n+1 is "air".  n=1 is the bottom of
     !  the atmosphere.
     ! Note:  media name of the node is obtained by
     !            Media(MediaNo) where MediaNo is fixed by
     !              MediaNo= atmos%node2mediaNo(node) just after
     !            the present active track is extracted by calling
     !     cpop.
     !            NoOfMedia (# of diff. emdia types in the current
     !            atmospher.  Medai, MediaNo, NoOfMedia are in
     !            ZmediaLoft.h
     real(8):: rhoc(maxnodes) ! rho correction factor give like H2O*0.01
     character(len=MAX_MEDIANAMELENG)::matter(maxnodes)
     
     integer:: nodes  ! actual # of total nodes.

  end type atmosph

  type(atmosph):: atmos

end module modAtmosDef

subroutine creadObjParam(io)
  use modAtmosDef
  implicit none
#include "Zmanagerp90.h"

  integer,intent(in):: io  ! logical dev # for ObjParam

!  real(8)::  EradiusIn  ! Earth radius beofre ObjParam read
  
  namelist /ObjParam/ &
   Eradius, &   ! radius of the Earth 6378136.0d0 m.
                ! sun--> 695510.0d3 m       
   Maxheight, & ! =210.d3 !  max height with some air. m 800d3 ? 
                !  3e6   sun
   Minheight, & ! -500.d0 !min height without soil. m
                ! -70.d0  sun
   MagChgDist, & ! =20.d3   !2 Distance where mag. can be seen
                 !  2000d3    sun ??
   gaccel,   &  !  = 9.80665d0   ! gravitational accel. m/sec2
                ! at near the surface of the Earth
                !  275.0    sun
   Hinf,   &  !  = 200.d3 ! m  At very high altitude
              !  = 2.d6  m  sun 
  !  air density is calulated using
  ! this value as the scale height.  rho= rho_m exp( -(z-Zm)/Hinf)
  !     -M means the value at Maxheight.
   AlmostVacH  !  = 210.0d3   ! m
               !  = 3e6  m  sun   
!   AlmostVacT   ! 1.8d-3  !  kg/m2  this will be fixed using AlmostVacH
  !
  
  if(io == 5 ) then
     read(*, ObjParam)
  else
     read(io, ObjParam)
  endif
  if( AtmosFile == " " ) then
     write(0,*) ' You gave non-blank "ObjFile"; it means you want to '
     write(0,*) ' define non Earth; in that case you have to give a value'
     write(0,*) ' of "AtmosFile" in "ObjParam" in "ObjFile" which is: '
     write(0,*)  trim(ObjFile)
     write(0,*) ' If you want to use Earth with some non standard atmospere'
     write(0,*) ' you may keep "OjbFile" blank and give a value to "AtmosFile"'
     stop
  endif
      
end subroutine creadObjParam

subroutine callocAtmosTbl
  use modAtmosDef
  implicit none
  
  
  allocate( ThickTbl(maxl2t))
  allocate( LenTbl(maxl2t))
  allocate( CosTbl(maxl2t))
  allocate( HeightTbl(maxl2t) )
end subroutine callocAtmosTbl
