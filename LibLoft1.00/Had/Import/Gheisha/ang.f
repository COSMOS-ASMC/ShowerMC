c
      subroutine ang(k,l,u,v)
c
c --- pv-array ---
      common /vecuty/ pv(10,200)
c
c
      double precision a,b,c,d
c
      a=0.0
      b=0.0
      c=0.0
      do 38 i=1,3
      a=a+pv(i,k)*pv(i,k)
      b=b+pv(i,l)*pv(i,l)
      c=c+pv(i,k)*pv(i,l)
 38   continue
      d=sqrt(a*b)
      if (d .ne. 0.0) d=c/d
      if (abs(d) .gt. 1.d0) d=sign(1.d0,d)
      u=d
      v=acos(d)
      return
      end
