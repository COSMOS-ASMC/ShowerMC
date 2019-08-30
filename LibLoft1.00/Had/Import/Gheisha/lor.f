c
      subroutine lor(k,l,m)
c
c --- pv-array ---
      common /vecuty/ pv(10,200)
c
c
      double precision a,b,c
c
      a=0.0
      do 6 i=1,3
      a=a+pv(i,k)*pv(i,l)
 6    continue
      a=(a/(pv(4,l)+pv(5,l))-pv(4,k))/pv(5,l)
      b=pv(5,k)*pv(5,k)
      do 8 i=1,3
      c=pv(i,k)+a*pv(i,l)
      b=b+c*c
      pv(i,m)=c
 8    continue
      pv(4,m)=sqrt(b)
      pv(5,m)=pv(5,k)
      return
      end
