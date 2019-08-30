subroutine kvec_prod3(v1, v2, v)
  real(8),intent(in):: v1(3), v2(3)
  real(8),intent(out):: v(3)

  v(1) = v1(2)*v2(3) - v1(3)*v2(2)
  v(2) = v1(3)*v2(1) - v1(1)*v2(3)
  v(3) = v1(1)*v2(2) - v1(2)*v2(1) 
end subroutine kvec_prod3


  


