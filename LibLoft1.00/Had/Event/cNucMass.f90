  module NucMassFormula
    implicit none
    real(8),parameter:: AtomicMassUnit = 931.494028d0 ! MeV
    integer,save:: formula=2
    contains
  !  give mass in MeV of a given A,Z ion
      function cNucMass(A,Z) result(mass)
        implicit none
        integer,intent(in):: A  ! mass #
        integer,intent(in):: Z  ! charge
    
        real(8):: mass  ! in  MeV.

        mass = cNucMassWOBE(A,Z) - cNucMassBE(A,Z)
      end function cNucMass

      function cNucMassWOBE(A,Z) result(mass)
        implicit none
        integer,intent(in):: A  ! mass #
        integer,intent(in):: Z  ! charge
        real(8):: mass    ! MeV
        mass =938.27231d0*Z + 939.56563d0*(A-Z)
      end function cNucMassWOBE

      function cNucMassBE(A,Z) result(BE)
        implicit none
        integer,intent(in):: A  ! mass #
        integer,intent(in):: Z  ! charge
        real(8):: BE    ! binding energy in MeV
        if( A <= 11 .and. Z<=5  ) then
           BE = cNucMassBEsmallA(A,Z)
        else
           BE = cNucMassBElargeA(A,Z)
        endif
      end function cNucMassBE

      function cNucMassBElargeA(A,Z) result(BE)
        implicit none
        integer,intent(in):: A  ! mass #
        integer,intent(in):: Z  ! charge
        real(8):: BE    ! binding energy in MeV
        if( formula == 1) then
           BE = cNucMassBE1(A,Z)
        elseif( formula == 2 ) then
           BE = cNucMassBE2(A,Z)
        else
           write(0,*) ' error formula= ', formula
        endif
      end function cNucMassBElargeA
      
      function cNucMassBEsmallA(A,Z) result(BE)
!          values after "else" is ????
        implicit none
        integer,intent(in):: A  ! mass #
        integer,intent(in):: Z  ! charge
        real(8):: BE    ! binding energy in MeV
        if(Z == 1 ) then
           if( A == 1 ) then
              BE = 0.
           elseif( A == 2 ) then
              BE = 1.112
           else
              BE = 1.2
           endif
        elseif(Z == 2 ) then
           if( A == 3 ) then
              BE = 2.57
           elseif( A == 4 ) then
              BE = 7.07
           else
              BE = 7.5
           endif
        elseif( Z == 3 ) then
           if( A == 6 ) then
              BE = 5.33
           elseif(A == 7) then
              BE = 5.61
           else
              BE = 6.0
           endif
        elseif( Z == 4 ) then
           if( A== 9 ) then
              BE = 6.46
           else
              BE = 6.9 
           endif
        elseif( Z == 5 ) then
           if( A == 10 ) then
              BE = 6.48
           elseif( A == 11 ) then
              BE = 6.93
           else
              BE = 6.93
           endif
        else
           write(0,*) ' A, Z=', A,Z
           write(0,*) ' invalid for cNucMassBEsmallA'
        endif
      end function cNucMassBEsmallA

      function cNucMassBE2(A,Z) result(BE)
!    gives binding energy, BE, in MeV.  shown in
!   J.K. Tuli, Nuclear Wallet Cards, National Nuclear Data
!   H.H. Heckman, PJ. Lindstrom, GH.D. Westfall and HJ. Center, Brook!   haven National Lab. (1985).
!  But the version is year 2005.
!  See
!     http://oregonstate.edu/instruct/ch374/ch418518/
        implicit none
        integer,intent(in):: A  ! mass #
        integer,intent(in):: Z  ! charge

        real(8):: BE    ! binding energy in MeV

        real(8):: delta   !
!      ci in MeV
        real(8),parameter::av=15.56, as=17.23, ac =0.7
        real(8),parameter::aa=23.285
        real(8):: N, A13
    
        N= A-Z
        A13 = dble(A)**0.33333333


        if( mod(Z,2) == 0 .and. mod(A,2) ==0 ) then
           delta =   11.0/sqrt(dble(A))    ! MeV
        elseif( mod(Z,2) /= 0 .and. mod(A,2) /= 0) then
           delta =  -11.0/sqrt(dble(A))   
        else
           delta = 0.
        endif

        BE=av*A -as*A13**2 -ac*Z**2/A13 - aa*(N-Z)**2/A + delta
      end function cNucMassBE2

      function cNucMassBE1(A,Z) result(BE)
!    gives binding energy, BE, in MeV.
!   By Myers and Swiatecki (Nucl. Phys. 81,1(1966).
!   as cited by 
!   J.K. Tuli, Nuclear Wallet Cards, National Nuclear Data
!   H.H. Heckman, PJ. Lindstrom, GH.D. Westfall and HJ. Center, Brook!   haven National Lab. (1985).
!  But the version is year 2005.
!  See
!     http://oregonstate.edu/instruct/ch374/ch418518/
        implicit none
        integer,intent(in):: A  ! mass #
        integer,intent(in):: Z  ! charge

        real(8):: BE    ! binding energy in MeV

        real(8):: asymm, delta   !
!      ci in MeV
        real(8),parameter::c1=15.667, c2=18.56, c3 =0.717
        real(8),parameter::c4=1.211, k=1.79
        real(8):: N, A13
    
        N= A-Z
        A13 = dble(A)**0.33333333
        asymm = 1-k*(dble(N-Z)/A)**2
        
        if( mod(Z,2) == 0 .and. mod(A,2) == 0 ) then
           delta =   11.0/sqrt(dble(A))    ! MeV
        elseif( mod(Z,2) /= 0 .and. mod(A, 2) /= 0 ) then
           delta =  -11.0/sqrt(dble(A))   
        else
           delta = 0.
        endif

        BE=c1*A*asymm -c2*A13**2*asymm -c3*Z**2/A13 + c4*Z**2/A + delta
      end function cNucMassBE1
      
      function cNucMassDefect(A,Z) result(deltaM)
        implicit none
        integer,intent(in):: A, Z
        
        real(8)::deltaM    !  in MeV
        
        deltaM = -cNucMassBE(A,Z)
      end function cNucMassDefect
        
      function cNucMassExcess(A,Z) result(delta)
!            mass - A in atomic mass unit
        implicit none
        integer,intent(in):: A, Z
        
        real(8)::delta    !  in MeV

        delta = cNucMass(A,Z)/AtomicMassUnit - A
      end function cNucMassExcess


      function cNucMassExcessMeV(A,Z) result(delta)
!            mass - A in MeV
        implicit none
        integer,intent(in):: A, Z
        
        real(8)::delta    !  in MeV

        delta = cNucMassExcess(A,Z)* AtomicMassUnit
      end function cNucMassExcessMeV

      function cProtonPhotoDissoci(A,Z) result(Gp)
        implicit none
        integer,intent(in):: A, Z
        real(8):: Gp
        Gp = min(dble(Z)/A, cgpByWH(Z))
      end function cProtonPhotoDissoci

      function cgpByWH(Z) result(Gp)
!        prob. of proton when photo-dissociation of charge Z
!        nuclei leads to the ejection of 1 nucleon.
!       E.V. Weinstock and J. Halpern, Phys. Rev. 94 (1954)
!      p.1651. 
!       Fitting by KK 
        implicit none
        integer,intent(in):: Z
        real(8):: Gp
        real(8):: x
        x = Z
        Gp = 0.967*exp(-(x/26.)**2.48)
      end function cgpByWH
      
      function cPhotoDissociWidth(A,Z) result(G)
!          only abundant ones
        implicit none
        integer,intent(in):: A, Z
        real(8):: G   !  MeV
        select case(A)
        case(12)    ! C
           G = 8.
        case(16)    ! O
           G = 10.
        case(18)    ! O
           G = 12.
        case(20:40) ! Ne, Si, S, Ar, Ca
           G = 10.
        case( 54 )  ! Fe
           G = 3.
        case( 56 )  ! Fe
           G = 5.
        case( 58 )  ! Ni 
           G = 10.
        case( 63 )  ! Cu
           G = 5.
        case( 90 )  ! Zr
           G = 4.0
        case( 107 ) ! Ag
           G = 5.0
        case( 160 ) ! Gd
           G = 4.0
        case( 197 ) ! Au
           G = 3.5
        case( 208 ) ! Pb
           G = 3.9
        case( 238 ) ! U
           G = 5.0
        case  default
           G = 10. 
        end select
      end function cPhotoDissociWidth
    end module NucMassFormula
