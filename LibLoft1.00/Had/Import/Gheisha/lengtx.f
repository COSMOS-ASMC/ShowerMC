c
      subroutine lengtx(k,u)
c
c --- pv-array ---
      common /vecuty/ pv(10,200)
c
c
      double precision a,b
c
      a=0.0
      do 36 i=1,3
      a=a+pv(i,k)*pv(i,k)
 36   continue
      b=sqrt(a)
      u=b
      return
      end
