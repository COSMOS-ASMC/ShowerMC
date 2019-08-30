#include "ZepMaxdef.h"
#include "ZepPos.h"

      integer maxvtx ! max number of vertex for one component
      real*8 gpsep   ! to show separater
      integer nvccl  ! a circle is approximated  by regular polygon of
                     ! nvccl vertexes.
      integer nvsph  ! to draw a sphere, polar angle is devided into
                     !   nvsph.
!      integer cut    !   0--> no cut in cyl, sphere
                     !   1--> cut along longitude in  cyl and sphere
                     !   2--> cut around the poles in  sphere
                     !   3--> cut around the poles and along 
                     !        longitude in sphere.
                     ! cyl means any cylindrical object such as
                     !     cyl ecyl, cone
                     ! sphere means  any spherical object such as
                     !     spherep, cap
      real*8 thetamin  ! works when cut=1 or 3. circle is not
      real*8 thetamax  ! drawn for azimuthal angle thetamin~thetamax
      integer maxcomp  ! max number of components
      integer*2 maxlevel ! max level given in config.
      integer maxdisplay ! max. num. of comps that can be displayed
                         !  without a lot of cpu time.
      integer*2 clevelmin   ! current min level of comp. to  be shown
      integer*2 clevelmax   ! current max level of comp. to  be shown
      integer*2 maxmaxlevel ! max level for which the number of comp. can be
                          !  counted.

      integer how      ! box surface display condition is
                       ! specified by 6 digit on/off.
      integer howcyl   ! 2 digits to specify drawing of floor and ceil of
                       ! cylinder like object.
      real*8 pamin     ! polar angle min (deg)
      real*8 pamax     ! polar angle max (deg)
      integer maxsel   ! selectable  nubmer of media/struc
      integer NsavedMX ! savable display cond.
      integer Nsaved   ! # of saved display conditions
      integer Ndflt    ! default display condition number
      integer Reverse  ! 0--> original, !=0--> Zrevmax-z --> z
      real*8  Zrevmax  ! max of config. used to reverse z
      integer confOrtrk ! 0-->config should be  first drawn
                        ! 1-->tracsk should be first drawn
      
      parameter ( maxsel = 100 )
      parameter ( maxmaxlevel = 8)
      parameter ( maxdisplay =MAX_DISPLAY )
      parameter ( gpsep = 2.7182d23 )
      parameter ( nvccl = 36 )
      parameter ( nvsph = 18 )
      parameter ( maxvtx = 50000 )
      parameter ( maxcomp = MAX_COMPONENT )
      parameter (NsavedMX = 5 )
      character*100 configFile
#if defined IBMAIX
      character*500 msg
#else
      character*1024 msg
#endif
      integer maxmat
      parameter (maxmat = 40)
      integer mcolor(maxmat)
      integer NoOfCompsEachMedia(maxmat)
      character*6 collist(maxmat)
      character*6 origcol(10)
      integer  maxspecifiable
      parameter (maxspecifiable = 100)
      integer*4 compnumb(maxspecifiable),
     *    mamordaught(maxspecifiable)
      character*8  matdef(maxmat)
      integer noOfMedia
      integer compnToMn(maxcomp)
      character*90 fn(maxmat)
      logical usedmedia(maxmat)
      integer offset
      parameter (offset=25)
      character*32  CondID(NsavedMX)
      character*32 CondFile
      integer  ahow(NsavedMX), ahowcyl(NsavedMX)
      real*8
     *  athetamin(NsavedMX), athetamax(NsavedMX),
     *  apamin(NsavedMX), apamax(NsavedMX)


      logical abyContain(NsavedMX), abyMedia(NsavedMX),
     *        abyStruc(NsavedMX), adrawall(NsavedMX)

      integer levelc(0:maxmaxlevel)  ! levelc(i) is the numb. of i-level comps.
      logical drawall, byContain, byMedia, byNumb,
     * byStruc
      integer drawmode, selmedia, containmode, selstruc
      integer maxmenu, nselstruc, nselMedia 
      integer subdnumber, NumComp,  NumWorld
      parameter (maxmenu=15)
      integer*4  draworhide(maxcomp)
      character*320 mlist, sellist
       type(epPos)::  pv(maxvtx)
      integer ndraw,  drawcomp(maxcomp)
      integer nstruc  ! no of defined structures in the config
      character*320 struclist, selstruclist
      character(MAX_STRUCCHR):: strucarray(maxsel), selstrucarray(maxsel), 
     * seledMedia(maxsel)
      common /zdrawc/  pv,
     *  athetamin, athetamax, apamin, apamax,
     * thetamin, thetamax, pamin, pamax, Zrevmax,
     * mcolor, subdnumber, NumComp, NumWorld,
     * how, howcyl, draworhide, noOfMedia, compnToMn, 
     * usedmedia, drawmode, drawall, adrawall, levelc,
     * selmedia, selstruc, containmode, ndraw,
     * nselstruc, nselMedia, nstruc, Nsaved,
     * drawcomp, byContain, byStruc, byMedia, byNumb,
     * compnumb, mamordaught, ahow, ahowcyl,
     * abyContain, abyMedia, abyStruc, Ndflt, Reverse,
     * confOrtrk, maxlevel, clevelmin, clevelmax,
     * NoOfCompsEachMedia
      common /zdawcc/CondID,  configFile, msg, mlist, sellist,
     *  matdef, fn, struclist, selstruclist, strucarray,
     *  selstrucarray, seledMedia, origcol, collist,
     *  CondFile



      

