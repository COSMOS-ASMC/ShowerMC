c
      subroutine sub(k,l,m)
c
c --- pv-array ---
      common /vecuty/ pv(10,200)
c
c
      double precision a,b
c
      a=pv(4,k)-pv(4,l)
      pv(4,m)=a
      b=a*a
      do 4 i=1,3
      a=pv(i,k)-pv(i,l)
      b=b-a*a
      pv(i,m)=a
 4    continue
      pv(5,m)=sign(sqrt(abs(b)),b)
      return
      end
