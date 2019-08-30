!      see max of  gzai(s) x G(s) or Psi(s) 
!          gzai(s) x G(s) is always < 1.0 
!          gzai(s) x Psi(s) = 1.003  or something like around s =0.5
!          so we may take that they are always <+ 1.0
!
       real*8 s
       real*8 gmigdl, psimig, gzai
       real*8 eps/1.d-4/, rho
       integer i
       rho = 1.e-3
       do i = 1, 5
          call zpart(7.25d0, 14.5d0, rho)
          do s =0.001, 1.0, 0.1
             write(*, *)sngl( s),sngl(  gmigdl(s, eps)* gzai(s))
             write(*, *)sngl( s) , sngl( psimig(s, eps)* gzai(s))
          enddo
          write(*,*)
          rho=rho/10.
       enddo
       end

