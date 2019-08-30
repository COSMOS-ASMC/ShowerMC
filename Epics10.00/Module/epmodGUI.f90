module modGUI
!      This is used for general userinerfaces which haven't
!   been forseen when the skeleton of this routine was
!   constructed.
  type gui
!       items for quenching effect follow
     real(8)::cf=1 !  output from epGUI. quenching factor to be multiplied to  qf*dE (dE for restricted energy loss)
     real(8)::qf=1 ! output from epGUI. quenching fraction
                   ! of dEdx (restricted energy loss rate)  
     real(8):: dedxf ! input
                 ! average dE/dx (GeV/(g/cm2)) for full energy loss
     real(8):: dedx  ! input
                 !  taverage dE/dx (GeV/(g/cm2)) for
                 ! restricted energy loss.
     integer:: modif ! input
           !if > 0, quenching is specified by modifier file. 
  end type gui

end module modGUI
