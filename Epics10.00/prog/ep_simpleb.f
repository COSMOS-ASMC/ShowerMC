!       ***************************
!          boundry for simple standard object. This is
!          formery in epbndry but this is accesssed by
!          say, angle and result in unresolved ref.
!          when drawConfig etc is used. so these are
!          separated.
        subroutine epbbox(comp, posli, dirli, el, jcon)
        implicit none
#include  "Zmedia.h"        
#include  "ZepTrackv.h"
#include  "Zcnfig.h"
!           find length to the boundary of comp component
!           from posl with direc cos dirl in local coord.
       type(Component):: comp
        integer jcon
        real*8 el, el2
!
        integer icon, icon2, icon3, icon4
       type(epPos)::   posli
       type(epDirec)::  dirli
        integer base
       type(epPos)::   posl
       type(epDirec)::  dirl
        
        base = comp%vol
        posl = posli
        dirl = dirli
        call kxplbx(
     *       posl%x, posl%y, posl%z, dirl%x, dirl%y, dirl%z,
     *       Volat( base+boxa), Volat(base+boxb), 
     *       Volat( base+boxc), el, jcon)
        return
!      ***********
       entry epbcyl(comp, posli, dirli,  el, jcon)
!      ***********
!      &&&&&&&&&&&&&&&&
       call epv2c_cyl(comp, posli, posl)
       call epv2cd_cyl(comp,dirli, dirl)
!       &&&&&&&&&&       
       call kxplcy(posl%x,  posl%y, posl%z, dirl%x,  dirl%y,
     *    dirl%z, Volat( comp%vol+cylr), Volat(comp%vol+cylh), 
     *    el, jcon, icon2)
        return
!      **********
       entry epbpip(comp, posli, dirli, el, jcon)
!

!      &&&&&&&&&&&&&&&&
       call epv2c_pipe(comp, posli, posl)
       call epv2cd_pipe(comp,dirli, dirl)
!       &&&&&&&&&&       

!             outer cylinder
       call kxplcy(posl%x, posl%y, posl%z, dirl%x, dirl%y, dirl%z, 
     *   Volat(comp%vol+pipeor), Volat( comp%vol+pipeh), 
     *    el, icon, icon2)
       if(icon2 .ne. -1)  then   !  there is crossing point
!                inner cylinder
          call kxplcy(posl%x, posl%y, posl%z, dirl%x, dirl%y, dirl%z, 
     *      Volat( comp%vol+pipeir), Volat(comp%vol+pipeh),
     *         el2, icon3, icon4)
          if(icon4 .ne.  -1) then
!                 there is crossing point with inner cyl. too
             el = min(el, el2)
             if(icon .eq. 0 .and. icon3 .eq. 1) then
                jcon = 0
             else
                jcon = 1
             endif
          else
!               xssing with outer cyl only.
             jcon = icon
          endif
       else
!             if there is no xssing point with outer cyl. no xsing at all
          jcon = -1
       endif
       return
!      **********
       entry epbsph(comp, posli, dirli, el, jcon)
!          compute  current pos to  comp sphere
!               current  pos is assumed to be inside the sphere
!     
       posl = posli
       dirl = dirli
       call kxplsph(posl%x, posl%y, posl%z,
     *      dirl%x, dirl%y, dirl%z, 
     *      Volat(comp%vol+sphr), el, jcon)
       return
       end
