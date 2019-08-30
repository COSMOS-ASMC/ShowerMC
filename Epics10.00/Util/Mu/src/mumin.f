!      mumin : minimum erergy to reach a given depth
      subroutine mumin(pathg, ec)
!       obtaine rough energy (ec) of muon which can run pathg (g/cm**2)
!          ec in tev.    muon of e<ec will run pathg with 1/2000 prob.
          x=log10(pathg)
          ec=10.**(  (x*0.7 - 6.25)*x + 13.2)
      end
