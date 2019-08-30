************************************************************************
*                                                                      *
      subroutine divcha(isc,line,nn)
*                                                                      *
*        devide succesive nonblank charactes one by one                *
*                                                                      *
*           isc :i  : scratch file unit number                         *
*           line:i/o: input data line                                  *
*           nn  :i  : length of input data line                        *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      character*(*) line

*-----------------------------------------------------------------------
*     re-try
*-----------------------------------------------------------------------

   10 continue

      if(line.eq.' ') goto 999

*-----------------------------------------------------------------------
*     find start point
*-----------------------------------------------------------------------

      ich1=0
      do 20 i=1,nn
         if(line(i:i).ne.' ') then
            ich1=i
            goto 30
         endif
   20 continue
      goto 999

*-----------------------------------------------------------------------
*     find end point
*-----------------------------------------------------------------------

   30 continue
      ich2=nn
      do 40 i=ich1,nn
         if(line(i:i).eq.' ') then
            ich2=i
            goto 50
         endif
   40 continue
   50 continue

*-----------------------------------------------------------------------
*     write successive character to scratch file
*-----------------------------------------------------------------------

      write(isc,'(a)') line(ich1:ich2)
      line(ich1:ich2)=' '
      goto 10
999   continue

*-----------------------------------------------------------------------

      return
       end subroutine


************************************************************************
*                                                                      *
      subroutine chnumb(hl,nn,icl1,icl2)
*                                                                      *
*     supply the starting and ending column number of non-blank        *
*     character.                                                       *
*     hl : searching character data                                    *
*     nn : total number of character data                              *
*     icl1 : starting column                                           *
*     icl2 : ending column                                             *
*     (created by k.kosako and modified by k.niita)                    *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      character hl*(*)

      character tub*1
      tub = char(9)

*-----------------------------------------------------------------------

      icl1 = 0
      icl2 = 0

      do 110 i = 1, nn

  110 if( hl(i:i) .ne. ' ' .and. hl(i:i) .ne. tub ) goto 120

      goto 900

  120 icl1 = i

      do 130 i = nn, icl1, -1

  130 if( hl(i:i) .ne. ' ' .and. hl(i:i) .ne. tub ) goto 140

  140 icl2 = i

  900 continue

      return
       end subroutine


************************************************************************
*                                                                      *
      function unirn(dummy)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
      unirn = fltrn(dummy)
      return
       end function


************************************************************************
*                                                                      *
      function fltrn(dummy)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common/irad/irands,irandf,nseed,ncall
      common /jcomon/ nabov,nobch,nocas,nomax
      common /rcomon/ rcasc
      data id1/1/
      data irandm/1/
c ----------------------------------------------------------------------
c
      if(irandm.eq.1) then
c        for mcnp-4a random generater .....

   10    fltrn=rang()
cKN      ncall=ncall+1

         if (fltrn.gt.1.0d+0) then
            write(6,20) fltrn,ncall,nocas,rcasc
            call parastop( 211 )
            go to 10
         endif

   20    format(' warning message from fltrn --- generated random ',
     &          'number was greater than 1.0, thus re-sample it.'
     &         /' rang=',e18.6,'  ncall=',i10,
     &          '   nocas=',i10,'  ncas=',f16.0)

      else

         if(id1.ne.1) go to 30
         if(irandf.eq.1) id1=irands
   30    id2=id1*663608941
         f2 = 0.5d+0 + id2 * 0.23283064d-9
         fltrn=f2
         f1=f2
         id1=id2
         nseed=id1
         ncall=ncall+1

      endif

      return
       end function



************************************************************************
*                                                                      *
      function gaurn(x)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
    1 y=exprnf(x)
      z=exprnf(x)
      if((y-1.)*(y-1.).gt.2.*z) go to 1
      if(unirn (dummy).lt.0.5) y=-y
       gaurn=y
      return
       end function


************************************************************************
*                                                                      *
      function exprnf(dummy)
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
      e=0.
    1 continue
      x= unirn(dummy)
      z=x
    5 y=unirn(dummy)
      if(z .le. y) go to 10
      z=unirn(dummy)
      if(z.lt.y) go to 5
      e=e+1.
      go to 1
   10 exprnf=e+x
      return
       end function


************************************************************************
*                                                                      *
      subroutine gtiso(u,v,w)
c ----------------------------------------------------------------------
c     calculate the unit directional vector of isotropic angular distr.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c ----------------------------------------------------------------------
c
    1 x=1.374736d+0*(unirn(dummy)-0.5d+0)
      y=1.374736d+0*(unirn(dummy)-0.5d+0)
      z=unirn(dummy)
      xsq=x*x
      ysq=y*y
      zsq=z*z
      d=xsq+ysq+zsq
      if(d*d.gt.z) go to 1
      u=2.d+0*x*z/d
      v=2.d+0*y*z/d
      w=(zsq-xsq-ysq)/d
      return
       end subroutine


************************************************************************
*                                                                      *
