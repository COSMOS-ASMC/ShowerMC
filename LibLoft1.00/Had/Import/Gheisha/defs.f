c
      subroutine defs(k,l,m)
c
c --- pv-array ---
      common /vecuty/ pv(10,200)
c
c
      double precision a,b
c
      mx=m
      my=m+1
      mz=m+2
      do 52 i=1,3
      f=pv(i,k)
      h=pv(i,l)
      pv(i,my)=f
      pv(i,mz)=h
 52   continue
      a=pv(2,my)*pv(3,mz)
      b=pv(3,my)*pv(2,mz)
      pv(1,mx)=a-b
      a=pv(3,my)*pv(1,mz)
      b=pv(1,my)*pv(3,mz)
      pv(2,mx)=a-b
      a=pv(1,my)*pv(2,mz)
      b=pv(2,my)*pv(1,mz)
      pv(3,mx)=a-b
      a=pv(2,mz)*pv(3,mx)
      b=pv(3,mz)*pv(2,mx)
      pv(1,my)=a-b
      a=pv(3,mz)*pv(1,mx)
      b=pv(1,mz)*pv(3,mx)
      pv(2,my)=a-b
      a=pv(1,mz)*pv(2,mx)
      b=pv(2,mz)*pv(1,mx)
      pv(3,my)=a-b
      do 58 j=mx,mz
      a=0.0
      do 54 i=1,3
      a=a+pv(i,j)*pv(i,j)
 54   continue
      b=sqrt(a)
      if (b .ne. 0.0) b=1.0/b
      do 56 i=1,3
      pv(i,j)=b*pv(i,j)
 56   continue
 58   continue
      return
      end
