!           life time in s
        real*8  t0pi, t0k, t0mu, t0ks, t0kl, t0pi0, t0d0, t0dc
        real*8  t0sigmap, t0sigmam, t0sigma0, t0gzai0, t0gzaim
        real*8  t0lambda, t0lambdac, t0eta, t0bomega
        real*8  t0ds, t0Xic, t0Xic0, t0omeC0,  t0tau, t0seethru
!           decay width in GeV
        real*8  wrho, womega, wphai, wDelta
!       change  2014.Sep. 
       parameter (
     1  t0pi=2.603d-8, t0k=1.238d-8, t0mu=2.197d-6,
!     1  t0pi=2.603d-8, t0k=1.2371d-8, t0mu=2.197d-6,
     2  t0ks=8.956d-11, t0kl=5.116d-8,       
!     2  t0ks=8.923d-11, t0kl=5.116d-8,           
     3  t0pi0=8.52d-17, t0d0=4.1d-13, t0dc=1.04d-12  )
!     3  t0pi0=8.4d-17, t0d0=4.2d-13, t0dc=1.06d-12  )
        parameter (
     1  t0sigmap = 0.8d-10,
     2  t0sigma0 = 5.8d-20,
     3  t0sigmam = 1.482d-10,
     4  t0gzai0 = 2.9d-10,
     5  t0gzaim = 1.641d-10,
     6  t0lambda = 2.623d-10,
     7  t0lambdac = 2.0d-13,
     8  t0eta     = 7.5d-19,
     9  t0bomega  = 0.822d-10 )

       parameter (
     4  wrho=150.d-3, womega=8.4d-3,
     5  wphai=4.4d-3, wDelta=112.5d-3    )
       real(8),parameter:: wetap=0.196d-3
       parameter (
     1  t0ds = 500.0d-15, t0Xic =442.0d-15,   t0Xic0 = 112.0d-15,
     2      t0omeC0 = 69.0d-15,  t0tau = 290.0d-15     )
       real(8),parameter:: t0etap=3.3d-21 ! =197.e-3*1.0e-13/weta/3e10)
       

