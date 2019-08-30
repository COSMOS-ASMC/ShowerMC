c
      subroutine trac(k,l,m)
c
c --- pv-array ---
      common /vecuty/ pv(10,200)
c
c
      double precision b,g(3)
c
      n=l
      do 62 j=1,3
      b=0.0
      do 60 i=1,3
      b=b+pv(i,n)*pv(i,k)
 60   continue
      g(j)=b
      n=n+1
 62   continue
      do 64 i=1,3
      pv(i,m)=g(i)
 64   continue
      return
      end
