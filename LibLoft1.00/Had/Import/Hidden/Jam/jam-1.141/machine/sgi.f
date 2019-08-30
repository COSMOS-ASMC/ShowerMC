************************************************************************
************************************************************************
*                                                                      *
*   PART 6: Machine dependent routines and utilities                   *
*                                           (For Silicon Graphic Cray) *
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      (S = subroutine, F = function, B = block data, E = entry)       *
*                                                                      *
*  S   JAMCPU   to measure and return the CPU time used                *
*  F   RN       to provide a random number generator                   *
*                                                                      *
************************************************************************
************************************************************************
      subroutine jamcpu(iunit,ic,iseed,isec)
************************************************************************
!  Purpose:   count the cpu time  U77 library
!  Variables: IC     -  0->start, 1->end
!  Function SECNDS returns the number of seconds that have elapsed
!    since midnight, less the value of its argument. The SECNDS routine
!    is useful for computing elapsed time of the execution of code.
!  The TIME subroutine returns the current system time.
!  For DTIME see  %man 3f dtime

      external dtime
      real tarray(2)
      character*24 fdate,today
      character timestr*8
      data cptime1/0.0/
      save cptime1, stime

      if (ic.eq.0) then
        cptime1=dtime(tarray)
        today = fdate()
        call time(timestr)
        write(iunit,'(''***********************************'')')
        write(iunit,'(''Starting time = '',a8,''  '',a9)') timestr,today
        stime=secnds(0.0)
        if (iseed.eq.0) then ! Set a different random seed each time 
          read(timestr,'(i2,1x,i2,1x,i2)') iseed1,iseed2,iseed3
          iseed = iseed1*10000 + iseed2*100 + iseed3
          if (iseed/2*2.eq.iseed) iseed = iseed + 1
        endif
        call srand(iseed) ! set a seed of randum number
        return

c...Compute Elapse and CUP time
      elseif (ic.eq.1) then
        cptime2=dtime(tarray)
        today = fdate()
        call time(timestr)
        write(iunit,'(''  ending time = '',a8,''  '',a9)')
     $  timestr, today
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
        return
      endif
      end

************************************************************************
      function rn(idumm)
************************************************************************
c...Purpose:  link random number generator.
      real*8 rn, rand
      rn = rand()
c     rn = pyrnd(0)    ! random number from pythia
      end
