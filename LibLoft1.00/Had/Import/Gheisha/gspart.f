      subroutine gspart(ipart,napart,itrtyp,amass,charge,tlife)
      implicit none
      integer *4 maxp
      parameter (maxp=200)
c.
c.    ******************************************************************
c.    *                                                                *
c.    *       store particle parameters                                *
c.    *                                                                *
c.    *    ==>called by : <user>, gpart                                *
c.    *       author    r.brun  *********                              *
c.    *                                                                *
c.    ******************************************************************
c.
      integer *4 ipart,itrtyp
      real *4 amass,charge,tlife
      character*(*) napart
      character *20 name
c
      integer *4 typp(maxp)
      real *4 amassp(maxp)
      real *4 chargep(maxp)
      real *4 tlifep(maxp)
      character *20 namep(maxp)
      save typp, amassp, chargep, tlifep, namep
      integer *4 fail
      character *10 c10
      integer *4 i,start,finish
      character*20 bl20
      save bl20
      data bl20/'                    '/
      data namep/maxp*'                    '/
c
      fail=1301
      if(ipart.le.0.or.ipart.gt.maxp) go to 1313
c
      typp(ipart)=itrtyp
      amassp(ipart)=amass
      chargep(ipart)=charge
      tlifep(ipart)=tlife
      namep(ipart)=napart
c
      return
c
      entry gfpart(ipart,name,itrtyp,amass,charge,tlife)
c
      fail=1302
      if(ipart.le.0.or.ipart.gt.maxp) go to 1313
c
       itrtyp=typp(ipart)
       amass=amassp(ipart)
       charge=chargep(ipart)
       tlife=tlifep(ipart)
       name=namep(ipart)
c
c      write(6,1991) ipart, name, amass, charge
c1991  format(1h , i3, a20, 3x, f10.5, 3x, f4.1) 
       return
c
      entry gpartp(ipart)
c
      fail=1303
      if(ipart.lt.0.or.ipart.gt.maxp) go to 1313
      start=1
      finish=maxp
      if(ipart.ge.1) then
       start=ipart
       finish=ipart
      endif !ipart.ge.1
       write(6,1002)
1002  format
     x(1x,' id name                    charge      mass  lifetime')
      do 1000 i=start,finish
       if(namep(i).eq.bl20) go to 1000
       if(tlifep(i).ge.1.0e+10) then
        c10='    stable'
       else
        write(c10,'(f10.3)') tlifep(i)
       endif !tlifep(i).ge.1.0e+10
       if(amassp(i).lt.1.0) then
        write(6,1003) i,namep(i),chargep(i),amassp(i),c10
       else
        write(6,1001) i,namep(i),chargep(i),amassp(i),c10
       endif !amassp(i).lt.1.0
1001   format(1x,i3,1x,a20,2f10.3,a10)
1003   format(1x,i3,1x,a20,f10.3,f10.6,a10)
1000  continue
c
      return
c
1313  continue
       write(6,*) 'fail=',fail,' in gspart'
       if(fail.eq.1301) then
        write(6,*) 'ipart=',ipart,' is out of range'
        write(6,*) 'this particle remains undefined'
       else if(fail.eq.1302) then
        write(6,*) 'ipart=',ipart,' is out of range'
        write(6,*) 'no valued returned of particle parameters'
       else if(fail.eq.1303) then
        write(6,*) 'ipart=',ipart,' is out of range'
        write(6,*) 'nothing will be printed'
       endif !fail.eq.1301
      return
      end
