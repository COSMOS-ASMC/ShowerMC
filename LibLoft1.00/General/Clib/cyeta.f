!         *********************************************
!           get rapidity and psude-rapidity from 
!           /ptcl/ data. 
          subroutine cyeta(p, y, eta)
!         *********************************************
          implicit none
!----          include '../Zptcl.h'
#include  "Zptcl.h"
          type(ptcl):: p
          real*8 y, eta
!
          real*8 pp, e, pt2, tm2, cost, sint, absp, tanh
           pp=p%fm%p(3)
           e=p%fm%p(4)
           pt2=p%fm%p(1)**2 + p%fm%p(2)**2
           tm2=pt2 + p%mass**2
           if(pp .ge. 0.d0) then
               y=log( (e+pp)**2/tm2 ) /2
           else
               y=log( tm2/(e-pp)**2 ) /2
           endif
           absp = sqrt(pt2 + pp**2)
           cost=pp/absp
           if(abs(cost) .gt. .99 ) then
                sint=sqrt(pt2) / absp
           else
                sint=sqrt(1.d0-cost**2)
           endif
           if(cost .gt. 0.d0) then
                tanh=sint/(1.d0+cost)
           else
                tanh=(1.d0 - cost)/sint
           endif
           if(tanh .le. 0.0) then
              if( cost .gt. 0.) then
                 eta = 20.
              else
                 eta = -20.
              endif
           else
              eta=-log(tanh)
           endif
       end







