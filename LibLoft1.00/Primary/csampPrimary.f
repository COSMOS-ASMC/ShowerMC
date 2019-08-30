!       *******************************************************
!       *  csampPrimary: samples a primary particle; its type
!       *                 and energy
!       *******************************************************
!   Note:  Here primary spectrum information is contained in common
!       variable Prim (done by init routine). 
!       If we  give it to csampPrimary0 we get
!       one sampled priamary of which energy etc is given in
!       variable in Prim.  Only /ptcl/ information is returned
!       to the caller.
!
        subroutine csampPrimary(p, fin)
!          p: /ptcl/  output.  energy, particle code, subcode,
!                              charge, mass are set.
!  

        implicit none
#include  "Zmanagerp.h"
#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryv.h"

        type(ptcl)::  p
        integer fin ! output. if no more sim. -->1 else 0


        fin = 0
        if( DestEventNo(2) .lt. 0) then
           if(abs(DestEventNo(2)) .le.  Prim%NoOfSamplings) then
              fin = 1
           endif
        endif
        if(fin .eq. 0) then
           call csampPrm0(Prim)


!          call cconv_prim_e(Prim)   ! to total energy
           call cconv_prim_e(Prim%each(Prim%label), Prim%sampled_e, 
     *      Prim%particle)         
           p = Prim%particle
           Prim%NoOfSamplings =  Prim%NoOfSamplings  + 1 ! update counter.
           Prim%NoOfSampComp(Prim%label, 1) =
     *          Prim%NoOfSampComp(Prim%label, 1) +1
        endif
        end
!       **********************   
        subroutine cupdtPrimC
        implicit none
#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryv.h"

        Prim%NoOfSampComp(Prim%label, 2) =
     *       Prim%NoOfSampComp(Prim%label, 2) +1
        end

!       ************************************
        subroutine csampPrm0(prm)
        implicit none

#include  "Zptcl.h"
#include  "Zprimary.h"
        type(primaries):: prm

        real*8 p_or_e
!
        integer i
        real*8 u
!
!          select one component
        call rndc(u)

        i = 1
        do while (u .gt. prm%cummInteFlux(i) )
           i = i + 1
        enddo
        prm%label = prm%each(i)%label
!
        call csampPrm1(prm%each(prm%label), p_or_e)


        prm%sampled_e = p_or_e
        end
        subroutine csampPrm1(each, p_or_e)
        implicit none

#include  "Zptcl.h"
#include  "Zprimary.h"
        type(component):: each
        real*8 p_or_e
!
        real*8 e1temp, ombeta, u
        integer i
!
        if( each%no_of_seg .eq. 0) then
           p_or_e = each%emin       ! = emax
        else
           call rndc(u)
           i = each%no_of_seg 
           do while (u .gt. each%norm_inte(i))
               i = i - 1
           enddo
!             use i-th segment
           if(each%diff_or_inte .eq. 'd') then    
               ombeta = (1.d0 - each%beta(i))           
               call rndc(u)
               if(abs(ombeta) .gt. 1.d-6) then
                  e1temp = each%energy(i)**ombeta
                  p_or_e =( u* (each%energy(i+1)**ombeta - e1temp) +
     *             e1temp )** (1.d0/ombeta)
               else    
                  p_or_e =( each%energy(i)/each%energy(i+1) )**u
     *                * each%energy(i+1)                  
               endif
           elseif(each%diff_or_inte .eq. 'i') then
               call rndc(u)
               p_or_e = each%energy(i)* (1.d0 - 
     *          u*(1.d0 - each%norm_inte(i+1)/each%norm_inte(i)))
     *          **(-1.D0/each%beta(i))
           else
               write(*, *) ' invlid diff_or_inte=',each%diff_or_inte,
     *         ' for primary=',each%symb
               stop 9999
           endif     
         endif  
         end
        subroutine cconv_prim_e(comp, e_or_p,  aPtcl)
        implicit none
!          given energy or momentum as given in primary 
!          specification is converted to total enregy.

#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zprimary.h"
         type(component):: comp  ! input. primary component
         real*8  e_or_p           ! input. energy or momentum as
                                  !        specified in primary data file
         type(ptcl)::aPtcl       ! output. 'partcle' of the primary 
         integer  code, massn
!       
         real*8 inGeV  ! p_or_e in  Gev
         real(8),parameter::hbarc=197.0e6*1e-15/1.e-9 ! hbarc in eV nm
         real(8),parameter::hc=hbarc*2*3.14159  ! hc/wl = eV
!
         inGeV= e_or_p *   comp%togev
         call cmkptc(comp%code,
     *               comp%subcode,
     *               comp%charge,
     *               aPtcl)        

         if(comp%code .eq. kmuon) then
!            give polarization for muon primary
            call csetMuonPol(0.d0)
         endif

         if(comp%etype .eq. 'e/p' .or.
     *      comp%etype .eq. 'e') then
            aPtcl%fm%p(4) = inGeV
     *                  
         elseif(comp%etype .eq.'p/p' .or.
     *      comp%etype .eq. 'p')then
            aPtcl%fm%p(4) = sqrt(aPtcl%mass**2 +
     *             inGeV**2)
         elseif(comp%etype .eq. 'ke/p' .or.
     *      comp%etype .eq. 'ke') then
            aPtcl%fm%p(4) = inGeV + aPtcl%mass
         elseif(comp%etype .eq. 'e/n') then
            code = comp%code
            if(code .eq. kgnuc) then
               aPtcl%fm%p(4) = aPtcl%subcode * inGeV
            elseif(code .ge. kalfa .and. 
     *         code .le. khvymax) then
               call cghvm(code, massn)
               aPtcl%fm%p(4) = massn * inGeV
            else
               aPtcl%fm%p(4) =  inGeV
            endif
         elseif(comp%etype .eq. 'ke/n') then
            code = comp%code
            if(code .eq. kgnuc) then
               aPtcl%fm%p(4) = aPtcl%subcode * inGeV +
     *            aPtcl%mass
            elseif(code .ge. kalfa .and. 
     *         code .le. khvymax) then
               call cghvm(code, massn)
               aPtcl%fm%p(4) = massn * inGeV +
     *            aPtcl%mass
            else
               aPtcl%fm%p(4) = inGeV +
     *            aPtcl%mass
            endif
	 elseif(comp%etype .eq. 'p/n') then
            code = comp%code
            if(code .eq. kgnuc) then
               aPtcl%fm%p(4) =
     *      	  sqrt( (aPtcl%subcode * inGeV) **2 + 
     *            aPtcl%mass **2)
            elseif(code .ge. kalfa .and. 
     *         code .le. khvymax) then
               call cghvm(code, massn)
               aPtcl%fm%p(4) =
     *      	  sqrt( (massn * inGeV) **2 + 
     *            aPtcl%mass **2)
            else
               aPtcl%fm%p(4) = sqrt(inGeV**2 +
     *            aPtcl%mass**2)
            endif
         elseif(comp%etype .eq. 'rig') then
            aPtcl%fm%p(4) =sqrt(
     *          ( inGeV*aPtcl%charge)**2 + aPtcl%mass**2 )
         elseif(comp%etype == 'nm' ) then
            if(comp%code /= klight ) then
               write(0,*) ' nm unit can be used only for light'
               stop
            endif
            aPtcl%fm%p(4) =  hc/e_or_p  ! hc/wl --> eV 
         else
            write(*, *) ' energy type=', comp%etype,
     *                  ' invalid. label=', comp%label, ' symb =',
     *                  comp%symb
            stop 9999
         endif
         if(comp%code == klight .and. comp%etype /= 'nm' ) then
            !  energy should be in eV so GeV to eV  conversion
            aPtcl%fm%p(4) = aPtcl%fm%p(4)*1.d9
         endif
        end
!       *******************************************************
        subroutine cconv_prim_e2(comp, E, e_or_p)
        implicit none
!         inverse of cconv_prim_e

#include  "Zptcl.h"
#include  "Zcode.h"
#include  "Zprimary.h"
         type(component):: comp  ! input. primary component
         real*8  E                ! input. total energy of a particle
                                  !   of 'comp' primary
         real*8 e_or_p            ! output. e or p (/n, etc) as given 
                                  !      in primary paticle file.
         integer  code, massn
!       
         real*8 Ex         !  e_or_p in GeV

         type(ptcl)::aPtcl

         call cmkptc(comp%code,
     *               comp%subcode,
     *               comp%charge,
     *               aPtcl)        


         if(comp%etype .eq. 'e/p' .or.
     *      comp%etype .eq. 'e') then
            Ex = E
         elseif(comp%etype .eq.'p/p' .or.
     *      comp%etype .eq. 'p')then
            Ex = sqrt(E**2 - aPtcl%mass**2)

         elseif(comp%etype .eq. 'ke/p' .or.
     *      comp%etype .eq. 'ke') then
            Ex = E - aPtcl%mass

         elseif(comp%etype .eq. 'e/n') then
            code = comp%code
            if(code .eq. kgnuc) then
               Ex = E/aPtcl%subcode
            elseif(code .ge. kalfa .and. 
     *         code .le. khvymax) then
               call cghvm(code, massn)
               Ex = E/massn
            else
               Ex = E
            endif
         elseif(comp%etype .eq. 'ke/n') then
            code = comp%code
            if( code .eq. kgnuc ) then
               Ex =( E -   aPtcl%mass)/aPtcl%subcode
            elseif(code .ge. kalfa .and. 
     *         code .le. khvymax) then
               call cghvm(code, massn)
               Ex =( E -   aPtcl%mass)/massn
            else
               Ex = E - aPtcl%mass
            endif
	 elseif(comp%etype .eq. 'p/n') then
            code = comp%code
            if( code .eq. kgnuc) then
               Ex =sqrt( E**2 -   aPtcl%mass **2)/aPtcl%subcode
            elseif(code .ge. kalfa .and. 
     *         code .le. khvymax) then
               call cghvm(code, massn)
               Ex =sqrt( E**2 -   aPtcl%mass **2)/massn
            else
               Ex = sqrt(E**2 - aPtcl%mass**2)
            endif
         elseif(comp%etype .eq. 'rig') then
            Ex = sqrt( E**2 - aPtcl%mass**2) /aPtcl%charge
         else
            write(*, *) ' energy type=', comp%etype,
     *                  ' invalid. label=', comp%label, ' symb =',
     *                  comp%symb
            stop 9999
         endif
        e_or_p = Ex/ comp%togev
        end
!       *******************************************************
!       *  cqPrimE: inquire sampled primary energy or p or rigidity
!       *                    as it is
!       *******************************************************
!
        subroutine cqPrimE(p_or_e)
!  
        implicit none

#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryv.h"
        real*8 p_or_e
        call cqPrimE0(Prim, p_or_e)
        end
!       *******************************************************
!       *  cqPrimE0: inquire sampled primary energy or p or rigidity
!       *                     as it is
!       *******************************************************
!
        subroutine cqPrimE0(prm, p_or_e)
!  
        implicit none

#include  "Zptcl.h"
#include  "Zprimary.h"
        real*8 p_or_e
        type(primaries):: prm
!
        p_or_e = prm%sampled_e
!
        end
!       *******************************************************
!       *  cqPrimLabel: inquire sampled primary label

!       *******************************************************
!
        subroutine cqPrimLabel(label)
!  
        implicit none

#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryv.h"
        integer label
        call cqPrimLabel0(Prim, label)
        end
!       ************************  inquire the number of primaries sampled
        subroutine cqNoOfPrim(no)
!       ************************
        implicit none
#include  "Zmanagerp.h"
        integer no   ! output.  no. of sampled primaries so far.
        no = EventNo
        end

!       *******************************************************
!
        subroutine cqPrimLabel0(prm,label)
!  
        implicit none
#include  "Zptcl.h"
#include  "Zprimary.h"
        integer label
        type(primaries):: prm

        label = prm%label
        end
!       *******************************************************
!           inquire all about current primaries
        subroutine cqPrimary(prm)
!         prm /primaires/ output.
!  
        implicit none
#include  "Zptcl.h"
#include  "Zprimary.h"
#include  "Zprimaryv.h"

        type(primaries):: prm
        prm = Prim
        end


