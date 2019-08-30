      subroutine cminTime2WebSec(obsdetxyz, ldep,  depidx,  awebmin)
      implicit none
#include  "Zglobalc.h"
#include  "Zmaxdef.h"
#include "Zmanagerp.h"
#include "Ztrack.h"
#include "Ztrackv.h"
#include "Zcode.h"
#include "Zheavyp.h"
#include "Zobs.h"
#include "Zobsp.h"
#include "Zobsv.h"
#include  "Zstackv.h"
#include "Zprivate.h"
      integer ldep ! layer number of the observation depth
      integer depidx ! layer index for web array
      type(coord):: obsdetxyz  ! input. observation detector axis in E-xyz
!                        = ObsSites(loc).pos.xyz  where  loc is the
!          observation layer number. 
      
      real*8 awebmin(nrbin, nfai, MaxNoOfSites) ! output. to store min time (ns)
                                 ! for each web sector.
      real*8 r, hdr, Ra, rc
      type(track):: inci
      type(coord):: angle
      integer i, j, ii, jj
      real*8 faimin,  fai, R0, Rbot, rmu, cosz
      data faimin/-15.0d0/
      type(coord):: xyz, oxyz, effpos
      real*8  temp,  dtemp, leng
      integer it



      call cqIncident(inci, angle)  ! we may better to use the first
                     ! collision point but not so easy to get it
                     ! since this is called before collision
      
      effpos = inci.pos.xyz    ! default is 100 km

!cc      effpos.z = 20.d3  ! so we use 20km instead
!         assume  web sector is aligned so that incident direcion
!         is on fai=0 of the web sector. 
!         sinz cosf= dirx
!         sinz sinf= diry
!         cosz = dirz 
!          tanf = diry/dirx;  fai=atan2(diry,dirx)
!          L dirx =x L diry =y  L dirz = z
!          L = z/dirz =
!      cosz = inci.vec.coszenith ~ -angle.r(3)

!                 Top view
!                       / this is web sector fai=0       
!                      / 
!                     /
!                    /      so fai* = fai-fai0
!                   / should be used as azimuthal angle of
!                  /           ptcls.  Therefore  we use
!                 /  fai0      effpos  as below
!               ------------------------> mag east
!
!      L cosz = effpos.z
!      L sinz = effpos.x
!      effpos.y = 0.
!

      leng = effpos.z/(-angle.r(3))
      effpos.x = sqrt(1.d0-angle.r(3)**2)*leng
      effpos.y = 0.
!          convert it to xyz system
      call cdet2xyz(obsdetxyz, effpos,  effpos)

      xyz.z = 0.
      xyz.x = 0.
      xyz.y = 0.
      call cdet2xyz(obsdetxyz, xyz, oxyz)
!/////////////
!      write(0,*)  ' layer =', ldep
!      write(0,*)  ' obsdetxyx=',obsdetxyz.x,
!     *       obsdetxyz.y,   obsdetxyz.z
!      write(0,*) ' center =', oxyz.x, oxyz.y, oxyz.z
!      write(0,*) ' effposx,y,x=',effpos.x,
!     *effpos.y, effpos.z
!//////////////

      R0 = sqrt( (oxyz.x-effpos.x)**2 +
     *           (oxyz.y-effpos.y)**2 +
     *           (oxyz.z-effpos.z)**2 )
!////////////
!      write(0,*) ' R0=', R0
!///////////

      call cgetMoliereU( ObsSites(ldep).pos.depth, cosz, rmu)
!//////////
!      write(0,*)  ' mu=', rmu
!///////////

      hdr = 10.**(bin/2.)

      do i = 1, nrbin
         rc = rbin(i) * rmu ! in m
         do j = 1, nfai
            Rbot = 10.d10
!               for all web sectors examine 4 corners
            do ii = 1, 2
!                rbin is the center of the web sector in r direction (log10 center).
              if(ii .eq.  1)  then
                 r = rc/hdr
              else
                 r = rc*hdr
              endif
              do jj = j, j+1
                 fai = faimin + (jj-1)*dfai
                 xyz.x = r*cos(fai*ToRad)
                 xyz.y = r*sin(fai*Torad)
                 call cdet2xyz(obsdetxyz, xyz, oxyz)
!                     *
!                    * ^
!                   *  ^
!                  *
!              R0 *     ^  Ra
!                *
!               *        
!              *         ^  
                 Ra = sqrt( (oxyz.x-effpos.x)**2 +
     *                 (oxyz.y-effpos.y)**2 +
     *                 (oxyz.z-effpos.z)**2 )
                  if( Ra-R0 .lt. Rbot ) Rbot = Ra-R0
               enddo
            enddo  ! for a given fai bi
!/////////////
!            if(i .gt. 20 .and. j. eq. 1) then
!               write(0,*) " ridx=",i, " Rbot=", Rbot
!            endif
!///////////
!/////////////
!            Rbot may not be real mininum if incident axis lies
!            on the sector. for safety  25 % correction.  in n sec
            temp = 0.35-0.1*(i/30.)
            dtemp =1.d9* ( Rbot - temp* abs(Rbot) )/c  ! ns

            awebmin(i,j,depidx) = dtemp   
!////////////////////
!            if(i .le. 3 .and. j .le. 2) then
!               write(0,*) 'i,j=',i,j, ' min t=', dtemp
!            endif
!///////////////
         enddo
      enddo


      end
