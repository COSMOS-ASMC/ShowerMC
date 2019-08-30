!               mass (GeV)

       real*8 &
       masele, maspic, maspi0, maskc, mask0, masd, &
       masmu, masp,  masn,  masrho, masomg, masphi, &
       maseta, masddb, masnnb, wrho,  womega, wphai, &
       massigmap, massigma0, massigmam, masgzai0,  &
       masgzaim, maslambda, maslambdac, masbomega, &
       masds, masXic, masXic0, masomC0, mastau

       
       parameter ( &
       masele=0.511d-3, maspic=139.5685d-3, maspi0=134.9642d-3, &
       maskc=493.646d-3, mask0=497.671d-3, masd=1869.d-3, &
       masmu=105.659d-3, masp=0.93827231,  masn=0.93956563, &
       masrho=768.d-3, masomg=782.d-3, masphi=1019.4d-3, &
       maseta=548.8d-3, masddb=2*masd, masnnb=2*masp, &
       wrho = 150.e-3, womega = 8.4e-3, wphai = 4.4e-3  )
       parameter(  &
      massigmap = 1.189, massigma0 = 1.192, massigmam=1.197, &
      masgzai0 = 1.314, masgzaim = 1.321,  &
      maslambda = 1.115, maslambdac = 2.282, masbomega=1.672) 
       parameter( &
       masds = 1.968, masXic = 2.468, masXic0 = 2.471, &
       masomC0 = 2.695,  mastau= 1.777  ) 
!        masddb and masnnb are the minimum value.
