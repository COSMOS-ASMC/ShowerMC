subroutine epcos2dir(cosa, cs, sn, w)
  ! from a given polar angle cos, make
  ! an azimuthally random direction cosine
  real(8),intent(in):: cosa
  real(8),intent(out):: cs, sn, w(3)
  ! cs, sn is uniform random azimuthal angle's  cos and sin
  !  w: direction cos vercor. w(3) may not be cosa if |cosa|>1.
  !  if so,  nearest 1 or -1 is employed.
  real(8):: sina,  cos

  cost = cosa
  sina = 1.d0 - cost**2

  if(sina < 0.d0 ) then
     sina = 0.
     cost =sign(1.0, cost)
  else
     sina = sqrt(sina)
  endif
  call kcossn(cs,sn)
  w(1) = sina * cs
  w(2) = sina * sn
  w(3) = cost
end subroutine epcos2dir
