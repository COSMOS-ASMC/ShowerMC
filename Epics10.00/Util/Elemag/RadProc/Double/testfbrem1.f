      common / bpcom /  al183z,e,ccz,emass,bcoef,fz ,z333
      write(*, *) ' enter Z, energy'
      read(*,*) z, ee
      e = ee
      dummy = zpart(z)
      v = 1.e-5
      vm = vmaxv(ee)
      do i =1, 1000
         f = fbrem(v)
         write(*,*) v, f, f*v
         v = v* 10.**0.02
         if(v .gt. vm) goto 10
      enddo
 10   continue
      end
