!c          to test the <v> sampled by m.c and expected value by formula
!c          fai(v)dv is the sampling function.
!c       <v> m.c = is from this distribution.
!c       <v> formula= int( v*fai(v) from epsilon to 1 )/ int(fai(v) from
!c                    epsilon to 1)
!c                    divident is int(0 to 1) - int(0 to eps)
!c                  and can be given be muparf, mubrmf, muhadf
!c                    devidor is the total cross-section
         ein= 10.
         epsb=0.01
         epsp=0.1
         epsh=0.1
!c          kgf rock
         z=12.8
         a=25.8
         zba=.493
         z2ba=6.30
         rho=3.02
!c          standard rock
         z=11.
         a=22.
         zba=.5
         z2ba=5.5
         rho=2.8
!c         frejus rock
         z=11.
         a=22.
         zba=.5
         z2ba=5.87
         rho=2.73
!c
       call mucset(z, a, zba, z2ba)
       call mupair(epsp)
       call mubrem(epsb)
       call muhadr(epsh)
       sum=0.
       n=100000
        do   i=1, n
          call muhadv(ein, v)
          sum=sum+v
        enddo
       vav=sum/n
       call muhadf(ein, v1)
       call muhadx(ein, xs1)
       call muhadr(1.)
       call muhadf(ein, v2)
       write(*,*) ' <v> by m.c.=',vav,' <v> by formula=',(v2-v1)/xs1
       end
