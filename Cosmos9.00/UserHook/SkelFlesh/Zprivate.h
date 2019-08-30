      integer NpMax
      parameter( NpMax = 60000 )
      integer  fnoB
      parameter  (fnoB=4)  ! hybrid AS output file no.

      structure /ob/
        integer*2 where, code, subcode, charge
        real atime, erg, mass, x, y, wx, wy, wz, zenith
      end structure

!      record /ob/o(NpMax)
	   type(ob)::o(NpMax)

      structure /parent/ 
         real*8 posx, posy, posz, coszenith
                                  ! these must be double
                                  ! otherwise small shift with
                                  ! children may appear
         real depth, colHeight
         real*8  height,  atime
         integer*2 where, code, asflag
         real  erg
      end structure
	type(parent)::p


      structure /child/
          integer*2 code, subcode, charge
          real fm(4), mass
      end structure

	type(child)::c

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



