
!               mass (GeV)

       real*8
     1  masele, maspic, maspi0, maskc, mask0, masd,
     2  masmu, masp,  masn,  masrho, masomg, masphi,
     3  maseta, masddb, masnnb, wrho,  womega, wphai,
     4  massigmap, massigma0, massigmam, masgzai0, 
     5  masgzaim, maslambda, maslambdac, masbomega,
     6  masds, masXic, masXic0, masomC0, mastau, masetap,
     7  masDelta
       
       parameter (
     1  masele=0.511d-3, maspic=139.5685d-3, maspi0=134.9642d-3,
     2  maskc=493.646d-3, mask0=497.671d-3, masd=1869.d-3,
     3  masmu=105.659d-3, masp=0.93827231,  masn=0.93956563,
     4  masrho=770.d-3, masomg=782.d-3, masphi=1019.4d-3,
     5  maseta=548.8d-3, masddb=2*masd, masnnb=2*masp,
     6  wrho = 150.e-3, womega = 8.4e-3, wphai = 4.4e-3  )
       parameter(
     1 massigmap = 1.189, massigma0 = 1.192, massigmam=1.197,
     2 masgzai0 = 1.314, masgzaim = 1.321, 
     3 maslambda = 1.115, maslambdac = 2.282, masbomega=1.672)
       parameter(
     1  masds = 1.968, masXic = 2.468, masXic0 = 2.471,
     2      masomC0 = 2.695,  mastau= 1.777,
     3      masetap =957.8e-3, masDelta=1.232
     4		 )
!        masddb and masnnb are the minimum value.
