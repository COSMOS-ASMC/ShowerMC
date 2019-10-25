!                Constants tables in primary class. 
!          PrimaryIdTbl is  set by calling csetPrimTbl
!          El is set by calling csetEUnitTbl
!
         type primaryid
           sequence
             character*16  symb              ! 'cno' etc
             integer code, subcode, charge      ! used to call cmkptc
         end type primaryid
         type eunit
           sequence
             character*3  symb               ! 'PeV' etc
             real*8       togev              !  energy * togev --> GeV energy
         end type eunit

         type(primaryid):: PrimaryIdTbl(NoOfSymbols)  ! all symbol list
         type(eunit):: ErgUnitTbl(maxErgUnit)  ! MeV etc to GeV
  
         common /Zprimaryc/ PrimaryIdTbl, ErgUnitTbl


