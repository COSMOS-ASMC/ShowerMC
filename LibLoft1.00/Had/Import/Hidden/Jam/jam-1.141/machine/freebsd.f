c***********************************************************************
c***********************************************************************
c                                                                      *
c        PART 6: Machine dependent routines for FreeBSD                *
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

c...Purpose: Count the cpu time  FreeBSD version.
c...Variables: IC - 0->start, 1->end

      character    today*26,chart(3)*17
      real*8 ar1(5),ar(5)
      dimension itoday(7)
      data ar1/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0/
      data chart/ 'total  cpu time =',
     $            'user   cpu time =',
     $            'system cpu time ='/
      save ar1

      if(ic.eq.0) then
        call dtime(ar1)
        call dateb(today,itoday)
        write(iunit,'(''******************************************'')')
        write(iunit,'(''starting time = '',a26)') today
c...Set a different random seed each time 
        if(iseed.eq.0) then
          iseed = itoday(1)*10000+itoday(2)*100+itoday(3)
          if(iseed/2*2.eq.iseed) iseed = iseed + 1
        end if
      else if(ic.eq.1) then

C...Compute Elapse and CUP time
        call dtime(ar)
        call dateb(today,itoday)
        write(iunit,'(''  ending time = '',a26)') today      
        stime=idtime()
        write(iunit,'(''  elapse time     = '',g15.4,'' sec'')') stime
        do i=1,3
          isec=ar(i)-ar1(i) 
          imin=isec/60
          ihrs=imin/60
          imin=imin-ihrs*60
          isec=isec-ihrs*3600-imin*60
          write(iunit,901)chart(i), ihrs,imin,isec
901       format(' * ',a17,i3,' h ',i3,' m ',i3,' s')
        end do
        isec=ar(1)-ar1(1) 
        write(iunit,'(''******************************************'')')
      else
        write(6,iunit)'(jamcpu:) ic should be 0 or 1 ic=',ic
        stop
      end if
 
      end
 
c***********************************************************************

      function rn(idumm)

c...Purpose: link random number generator.
      real*8 rn,pyrnd
      common/rseed/iseed
      save /rseed/

c...From Pythia routine.
      rn = pyrnd(0)
      end

