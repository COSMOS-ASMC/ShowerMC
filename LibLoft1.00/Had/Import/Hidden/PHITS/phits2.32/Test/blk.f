      block data blkdt
      real(8)::x, y
      common /abc/ x, y(10)
      data y(1)/5.0/
      data y(2)/-5.5/
      end block data blkdt



