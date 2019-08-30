c...User defined analyze routine.
c***********************************************************************
c...This routine will be called in every collision and decay.
      subroutine jamanaus(indd,npart)

c...Analysize collision spectra.
      include 'jam1.inc'
      include 'jam2.inc'
      dimension indd(100)

c...Useful vectors for analysis.
      ichanel=mste(1)
      icltyp=mste(2)

      kf1=kcp(1,1)   ! PDG code for colliding particle 1
      kf2=kcp(1,2)   ! PDG code for colliding particle 2
      kf3=kcp(2,1)   ! PDG code for outgoing  particle 1
      kf4=kcp(2,2)   ! PDG code for outgoing  particle 2
      i1=mste(21)    ! line number of the ingoing particle 1
      i2=mste(23)    ! line number of the ingoing particle 2
      i3=mste(25)    ! line number of outgoing particle 1
      i4=mste(27)    ! line number of outgoing particle 2

      end
