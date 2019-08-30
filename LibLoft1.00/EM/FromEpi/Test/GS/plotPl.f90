program main
implicit none
real(4)::kpnx
real(4)::cost
real(8)::teta
integer::i, l
real(8),parameter:: pi=asin(1.d0)*2


do l = 0, 10
   teta =0.
   do while (teta< 2*pi)
      cost=cos(teta) 
      write(*,'(i3, 1p,3g13.4)') l, teta, kpnx(l,cost), cost
      teta = teta + 0.05d0
   enddo
enddo
end program main



