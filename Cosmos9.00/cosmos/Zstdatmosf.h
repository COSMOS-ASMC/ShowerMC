!        to be included just before the execution code
!      density as a function of height.
      real* 8 fd1, fd0      
      real*8 a, z0, h0, zz
!            scale height = linear func of height, zz
      fd1(zz, a, z0,  h0 ) = (1.d0 + a* (zz - z0)/h0) **(-1.d0/a)
!            const scale height
      fd0(zz, z0, h0) = exp(- (zz -z0)/h0)
!--------------------
