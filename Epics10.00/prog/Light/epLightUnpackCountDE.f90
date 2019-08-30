!   To decompose CountDE = singed Bmnd. 
!   after extraction only d inherit negative sign
!    epgetmn, epgetd, epgetB: for getting mn, d, B
!    epLightUnpackCountDE: to get mn, d, B at once
integer function epLightGetmn(countde)
  implicit none
!        extract mn part from sBmnd in countde                                
  integer(2)::countde  ! input. CountDE                                      
  integer::temp1, temp2
                              !  countde  
                              !  3112   -3512
  temp1 = abs(countde)/10     !   311     351
  temp2 = temp1/100           !     3       3
  epLightGetmn = temp1 - temp2*100 !    11      51
end function epLightGetmn

integer function epLightGetd(countde)
  implicit none
!        extract d part from sBmnd in countde                                
  integer(2)::countde  ! input. CountDE                                      
  integer::temp1, temp2, icountde 
                              !  countde  
  temp2 = abs(countde)             !  3112   -3512
  temp1 = temp2/10     !   311     351
  temp1 = temp2 - temp1*10    !     2       2
  !  epLightGetd = sign(temp1, countde) !   2   -2 NG for Gfort
  icountde = countde
  epLightGetd = sign(temp1, icountde) !    
end function epLightGetd

integer function epLightGetB(countde)
  implicit none
!        extract B part from sBmnd in countde                                
  integer(2)::countde  ! input. CountDE                                      
  integer::temp1, temp2
                              !  countde  
                              !  3112   -3512
  temp1 = abs(countde)/1000    !    3        3
  epLightGetB = temp1
end function epLightGetB

subroutine epLightUnpackCountDE(countde, d, mn, B)
  implicit none
  integer(2)::countde  ! input CountDE:  signed Bmnd
  integer::d, mn, B    ! ouput

  integer temp1, temp2, temp3, icountde
                       ! countde   
                       !  3102    -3512
  temp1 = abs(countde)     !  3102     3512
  temp2 = temp1/10         !   310      351
  temp3 = temp2/100        !     3        3

  !  d =sign(temp1- temp2*10, countde)  ! NG for Gfort
  icountde = countde
  d =sign(temp1- temp2*10, icountde)
  mn = temp2 - temp3*100
  B = temp3
end subroutine epLightUnpackCountDE



