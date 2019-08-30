!     cepos2cos:    epos code to cosmos code conversion
!     ccos2epos:    cosmos code to epos code conversion
subroutine cepos2cos(id, code, subcode, charge)
  implicit none
  integer,intent(in):: id  ! epos particle  code
  integer,intent(out):: code ! cosmos code
  integer,intent(out):: subcode ! cosmos subcode
  integer,intent(out):: charge ! cosmos charge
      
  integer:: idtrafo  ! epos<--->pdg in epos-ids-199.f 
      
  integer:: idpdg  !  pdg code

  idpdg = idtrafo("nxs", "pdg", id)  ! id of epos is converted into pdg code
  
  call ckf2cos(idpdg, code, subcode, charge)
end subroutine cepos2cos

subroutine ccos2epos(code, subcode, charge, id)
  implicit none
  integer,intent(in):: code ! cosmos code
  integer,intent(in):: subcode ! cosmos subcode
  integer,intent(in):: charge ! cosmos charge
  integer,intent(out):: id  ! epos particle  code      

  integer:: idtrafo  ! epos<--->pdg in epos-ids-199.f 
      
  integer:: idpdg  !  pdg code
  
  call ccos2kf(code, subcode, charge, idpdg)
  id = idtrafo("pdg", "nxs", idpdg)  ! idpdf of pdg is converted into epos code
!////////////
!  if( code >= 9 ) then
!     write(0,'(a,i5, a, i4, a, i3)') &
!       ' code=',code, ' subc=',subcode, ' charege=',charge
!     write(0,*) ' are converted to pdgcode =', idpdg
!     write(0,*) ' is converted to epos code =', id
!  endif
!/////////////////
end subroutine ccos2epos

