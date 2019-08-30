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
c...In HP-ux,
c...Program must be compiled with either +e or +E1 option
c...in order to use  DATE, TIME,SECNDS and
c...+U77 option (BSD 3F library libU77) for DTIME
c...call date_and_time(datestr)
c...character*9 datestr = dd-mmm-yy
c...function SECNDS returns the number of seconds that have elapsed
c...since midnight, less the value of its argument. The SECNDS routine
c...is useful for computing elapsed time of the execution of code.
c...The TIME subroutine returns the current system time.
c...For DTIME see  %man 3f dtime

      external dtime
      real tarray(2)
      dimension idati(6)
      character today*9,timestr*8
      data cptime1/0.0/
      save  stime,cptime1

      if(ic.eq.0) then
        cptime1=dtime(tarray)
        call date_and_time(today)
        call pytime(idati)
        call time(timestr)
        write(iunit,'(''***********************************'')')
        write(iunit,'(''Starting time = '',A8,''  '',A9)')timestr,today
        stime=secnds(0.0)
c...Set a different random seed each time 
        if(iseed.eq.0) then
          read(timestr,100) iseed1,iseed2,iseed3
  100     format(i2,1x,i2,1x,i2)
          iseed = iseed1*10000 + iseed2*100 + iseed3
          if(iseed/2*2.eq.iseed) iseed = iseed + 1
        end if
      end if

C...Compute Elapse and CUP time
      if(ic.eq.1) then
        cptime2=dtime(tarray)
c       call date_and_time(today)
        call time(timestr)
        write(iunit,'(''  ending time = '',a8,''  '',5(i2,1x))')timestr
     $    ,idati(3),idati(1),idati(4),idati(5),idati(6)

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
        write(iunit,'(''***********************************'')')
        isec=cptime2-cptime1
      end if

      end

c***********************************************************************

      function rn(idumm)

c...Purpose:  link random number generator.
      real*8 rn,pyrnd
      common/rseed/iseed
      save /rseed/
c...Program must be compiled with either +e or +E1 option
c...in order to use  function RAN on HP.
c...for HP
c     rn = ran(iseed)
      rn = pyrnd(0)    ! random number from pythia
      end
