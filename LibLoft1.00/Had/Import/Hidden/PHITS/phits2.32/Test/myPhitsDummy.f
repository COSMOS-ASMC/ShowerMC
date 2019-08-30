      subroutine myPhitsDummy
      implicit real*8 (a-h,o-z)
#include "bert.inc"
#include "jam2.inc"
#include "param00.inc"
#include "param01.inc"
#include "param02.inc"
c         bert-bl0
      common/rtcom/tapcrs(6600)
      common/rtco1/tapcr1(352),tapcr2(352)
c         bert-bl1
      common /geosig/ geosig(250)
c          param00
      common /natnuc/ natnn(maxpt), natnm(maxpt,10), patnn(maxpt,10)
      common /maspn0/ nzz0(0:nnuc), nnn0(0:nnuc), jgs0(0:nnuc)
c           param01
       common /input1/ mstq1(mxpa1), parq1(mxpa1)
c         bert-bl0
      write(0,*) 'tapcrs '
      write(0,*) tapcrs(1:10)
      write(0,*) 'tapcr1 '
      write(0,*) tapcr1(1:10)
      write(0,*) 'tapcr2 '
      write(0,*) tapcr2(1:10)

c         bert-bl1
      write(0,*) 'geosig '
      write(0,*) geosig(1:10)
c         bert-bl2
      write(0,*) 'bert-bl2: sf etc '
      write(0,*)  sf, nrt, ln, part,  pnms,  dncms,  sqnm, rcpmv,
     * poms
      write(0,*) 'crsc'
      write(0,*)  crsc(1:10)
c         jam2.inc
      write(0,*) ' jam2: kchg'
      write(0,*)  kchg(1:10,1)
      write(0,*)  kchg(1:10,3)
      write(0,*)  kchg(1:10,5)

      write(0,*) 'natnn '
      write(0,*) natnn(1:10)
      write(0,*) ' natnm '
      write(0,*) natnm(1:10,1)
      write(0,*) ' patnn '
      write(0,*) patnn(1:10,1)

      write(0,*) 'nzz0'
      write(0,*) nzz0(0:9)

      write(0,*) ' mstq1, parq1 '
      write(0,*) mstq1(1:10)
      write(0,*) parq1(1:10)
      end




