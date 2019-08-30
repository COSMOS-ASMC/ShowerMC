      subroutine mudedx(e, dedt)
!          get muon energy loss rate tev / (g/cm**2) due to the
!         ionization, knockon.
!
       include  'Zmucom.f'
!           0.0109= mu**2/emass/2 in tev
        emu= e**2/(e+0.010924)
        dedt=zba*(3.776e-6 + 0.1535e-6* log(emu/mu))
      end
