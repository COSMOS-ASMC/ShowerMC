      integer NpMax
      parameter( NpMax = 200000 )
      integer  fnoB
      parameter  (fnoB=4)  ! hybrid AS output file no.
      type ob
        integer*2 where, code, subcode, charge
	real*8 user
        real atime, erg, mass, x, y, wx, wy, wz, zenith
      end type ob
	character*8 proc
      type(ob)::o(NpMax)

      type parent 
         real*8 posx, posy, posz, coszenith
                                  ! these must be double
                                  ! otherwise small shift with
                                  ! children may appear
         real depth, colHeight
         real*8  height,  atime
!///////////
	 real*8  user
!////////////
         integer*2 where, code, asflag
         real  erg
      end type parent
      type(parent):: p


      type child
          real*8 user
          integer*2 code, subcode, charge
          real fm(4), mass
      end type child

      type(child):: c

      integer Np, Wdev, Mdev, NhMin, NgMin, NoOfLowE, NLowCounter,
     *        Rdev, HowFlesh
      real*8 SumehMin, SumegMin,  Ethresh, Cuteg, Cutneg

      character*120  Mskel, Wskel

      integer Accepted, Where
      logical TopOfNode, RealBegin, RealEnd, Copy
!      save o, c, p, TopOfNode, RealBegin, RealEnd, Copy
!       save  Accepted, Where
      common /skelcom/ o, c, p,   SumehMin, SumegMin, Ethresh,
     *     Cuteg, Cutneg,
     *     Np, NoOfLowE, NLowCounter,
     *     Wdev, Mdev,  Rdev,
     *     NhMin, NgMin, Accepted, Where,  HowFlesh, TopOfNode,
     *     RealBegin, RealEnd, Copy

      common /cskelcom/ Wskel, Mskel

      logical bgevent
      type(ptcl):: cmsptcl,  incip, pjcm
      type(track):: incident     
      common /cother/ proc, incident, cmsptcl, incip, pjcm, bgevent
