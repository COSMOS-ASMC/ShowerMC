c***********************************************************************
c***********************************************************************
c                                                                      *
c        PART 6: Machine dependent routines and utilities (U77library) *
c                                                                      *
c   List of subprograms in rough order of relevance with main purpose  *
c      (S = subroutine, F = function, B = block data, E = entry)       *
c                                                                      *
c  s   jamcpu   to measure and return the CPU time used                *
c  f   rn       to provide a random number generator                   *
c                                                                      *
c***********************************************************************
c***********************************************************************
c***********************************************************************

      subroutine jamcpu(iunit,ic,iseed,isec)

c...Purpose:   count the cpu time  U77 library
c...Variables: IC     -  0->start, 1->end
c...call date_and_time(datestr)
c...function SECNDS returns the number of seconds that have elapsed
c...since midnight, less the value of its argument. The SECNDS routine
c...is useful for computing elapsed time of the execution of code.

c     external dtime
      real tarray(2)
      character today*8,gtime*10,gzone*5
      integer value(8)
      real tresult
      character timestr*8
      data cptime1/0.0/
      data timestr/'today  '/
      save  stime,cptime1

      if(ic.eq.0) then
        call dtime(tarray,tresult)
        cptime1=tarray(1)
        call date_and_time(today,gtime,gzone,value)
        write(iunit,'(48(''*''),)')
        write(iunit,800)value(1),value(2),value(3),value(5),value(6),
     $           value(7),gzone
800     format(' Starting time = ',i4'/',i2,'/',i2,'   ',
     $         i2':',i2,':',i2,'s   ',A5)

        stime=secnds(0.0)

c...Set a different random seed each time 
        if(iseed.eq.0) then
          read(gtime,100) iseed1,iseed2,iseed3
  100     format(i2,2x,i2,2x,i2)
          iseed = iseed1*10000 + iseed2*100 + iseed3
          if(iseed/2*2.eq.iseed) iseed = iseed + 1
        end if
      end if

C...Compute Elapse and CUP time
      if(ic.eq.1) then
        call dtime(tarray,tresult)
        cptime2=tarray(1)
        call date_and_time(today,gtime,gzone,value)
        write(iunit,801)value(1),value(2),value(3),value(5),value(6),
     $           value(7),gzone
801     format('   ending time = ',i4'/',i2,'/',i2,'   ',
     $         i2':',i2,':',i2,'s   ',A5)


        stime=secnds(stime)
        isec=stime
        imin=isec/60
        ihrs=imin/60
        imin=imin-ihrs*60
        isec=isec-ihrs*3600-imin*60
        write(iunit,901)ihrs,imin,isec
901     format(' * Elapse time =',i3,' h ',i3,' m ',i3,' s')

        isec=cptime2-cptime1
        imin=isec/60
        ihrs=imin/60
        imin=imin-ihrs*60
        isec=isec-ihrs*3600-imin*60
        write(iunit,902)ihrs,imin,isec
902     format(' *    CUP time =',i3,' h ',i3,' m ',i3,' s')
        write(iunit,'(48(''*''),)')
        isec=cptime2-cptime1
      end if

      end

c***********************************************************************

      function rn(idumm)

c...Purpose:  link random number generator.
      real*8 rn,pjrnd
      common/rseed/iseed
      save /rseed/
c...Program must be compiled with either +e or +E1 option
c...in order to use  function RAN on HP.
c...for HP
c     rn = ran(iseed)
      rn = pjrnd(0)    ! random number from pythia
      end
