!           common parameters for tracking
!
!	(->  -------------------------------------

        real*8  Deltpp  !2  p-p xsection increases as $E^{Deltpp}$(E$>$ 100GeV)
        real*8  Deltpip !2  pi-p xsection increases as $E^{Deltpip}$ (E$>$ 100GeV)
        real*8  Deltkp  !2  k-p xsection increases as $E^{Deltkp}$ (E$>$ 100GeV)
        real*8  IncreaseXsec !2   how the xsection increases. 1.0$\rightarrow$  power of E
!              above ones are obsolete
!	<-)  -----------------------------------

        common /Zxsectionp/
     *  Deltpp, Deltpip, Deltkp, IncreaseXsec
!

