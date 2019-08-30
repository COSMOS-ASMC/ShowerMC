      FUNCTION BSPFIT(N,CF,X0,DX,XLOG,YLOG,XX)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
! Give function value and partial derivatives
!
!            ix+1
!            ___ 
!  F(X,Y) =  >    c(i)*Bspl(i,x1)
!            --- 
!           i=ix-2
!
!    Where ix=int((x-x0)/dx),
!          x1= (x-x0)/dx + 1
!
! Fitting parameters
      DIMENSION CF(N)
      LOGICAL XLOG,YLOG, newcoef
! Input
!    x: coordinate of the point
! Working area
      DIMENSION C0(4,4), CC0(4)
      data c0/-1.d0, 3.d0,-3.d0, 1.d0, 3.d0,-6.d0, 0.d0, 4.d0,
     &        -3.d0, 3.d0, 3.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
      save c0, cf1, cc0, ix0
      data ix0/-1000000000/
       IF (XLOG) THEN 
         X=LOG10(XX)
      ELSE
         X=XX
      END IF
       X1 = (X-X0)/DX + 1.0D0
      IX = INT(X1+100000.0) -100000
       newcoef = (ix.ne.ix0).or.(cf1.ne.cf(1))
       IBGN = MAX(1,IX-1)
      IEND = MIN(N,IX+2)
       s = x1 - real(ix)
      if(newcoef) then
         do k=4,1,-1
            cc0(k) = 0.
            DO 3 I=IBGN, IEND
               J = I - IX + 2
               cc0(k) = cc0(k) + cf(i)*c0(k,j)
 3          CONTINUE
         enddo
         cf1 = cf(1)
         ix0 = ix
      end if
       b = cc0(4) + s*(cc0(3)+s*(cc0(2)+s*cc0(1)))
       IF (YLOG) THEN
         BSPFIT  = EXP(B*LOG(10.D0))
      ELSE
         BSPFIT  = B
      END IF
!
      RETURN
      END
!
      FUNCTION BSFITD(N,CF,X0,DX,XLOG,YLOG,XX)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
! Give function value and partial derivatives
!
!            ix+1
!            ___ 
!  F(X,Y) =  >    c(i)*Bspl(i,x1)
!            --- 
!           i=ix-2
!
!    Where ix=int((x-x0)/dx),
!          x1= (x-x0)/dx + 1
!
! Fitting parameters
      DIMENSION CF(N)
      LOGICAL XLOG,YLOG,newcoef
! Input
!    x: coordinate of the point
!     Working area
      dimension cc0(4), cc1(3)
      dimension c0(4,4), c1(3,4)
      data c0/-1.d0, 3.d0,-3.d0, 1.d0, 3.d0,-6.d0, 0.d0, 4.d0,
     &        -3.d0, 3.d0, 3.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
       data c1/-3.d0, 6.d0,-3.d0, 9.d0,-12.d0, 0.d0, 
     &        -9.d0, 6.d0, 3.d0, 3.d0, 0.d0, 0.d0/
      save c0, c1, cf1, cc0, cc1, ix0
      data ix0/-1000000000/
       IF (XLOG) THEN 
         X=LOG10(XX)
      ELSE
         X=XX
      END IF
      X1 = (X-X0)/DX + 1.0D0
      IX = INT(X1+100000.0) -100000
       newcoef = (ix.ne.ix0).or.(cf1.ne.cf(1))
       IBGN = MAX(1,IX-1)
      IEND = MIN(N,IX+2)
       s = x1 - real(ix)
      if(newcoef) then
         do k=1,3
            cc1(k) = 0.
            DO I=IBGN, IEND
               J = I - IX + 2
               cc1(k) = cc1(k) + cf(i)*c1(k,j)
            enddo
         enddo
         cf1 = cf(1)
         ix0 = ix
      end if
      db = (cc1(3)+s*(cc1(2)+s*cc1(1)))
       BSFITD  = DB/DX
      IF (YLOG) THEN
         if(newcoef) then
            do k=4,1,-1
               cc0(k) = 0.
               DO I=IBGN, IEND
                  J = I - IX + 2
                  cc0(k) = cc0(k) + cf(i)*c0(k,j)
               enddo
            enddo
         end if
          b = cc0(4) + s*(cc0(3)+s*(cc0(2)+s*cc0(1)))
         BSFITD  = LOG(10.D0)*BSFITD*EXP(B*LOG(10.D0))
      END IF
      IF (XLOG) THEN
         BSFITD  = BSFITD/XX/LOG(10.D0)
      END IF
      RETURN
      END


      subroutine BsFITvdv(N,CF,X0,DX,XLOG,YLOG,XX, v, dv)
      IMPLICIT REAL*8 (A-H,O-Z)
! Give function value and partial derivatives
!
!            ix+1
!            ___ 
!  F(X,Y) =  >    c(i)*Bspl(i,x1)
!            --- 
!           i=ix-2
!
!    Where ix=int((x-x0)/dx),
!          x1= (x-x0)/dx + 1
!
! Fitting parameters
      DIMENSION CF(N)
      LOGICAL XLOG,YLOG,newcoef
! Input
!    x: coordinate of the point

      dimension cc0(4), cc1(3)
      dimension c0(4,4), c1(3,4)
      data c0/-1.d0, 3.d0,-3.d0, 1.d0, 3.d0,-6.d0, 0.d0, 4.d0,
     &        -3.d0, 3.d0, 3.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/

      data c1/-3.d0, 6.d0,-3.d0, 9.d0,-12.d0, 0.d0, 
     &        -9.d0, 6.d0, 3.d0, 3.d0, 0.d0, 0.d0/
      save c0, c1, cf1, cc0, cc1, ix0
      data ix0/-1000000000/

      IF (XLOG) THEN 
         X=LOG10(XX)
      ELSE
         X=XX
      END IF
      X1 = (X-X0)/DX + 1.0D0
      IX = INT(X1+100000.0) -100000

      newcoef = (ix.ne.ix0).or.(cf1.ne.cf(1))

      IBGN = MAX(1,IX-1)
      IEND = MIN(N,IX+2)

      s = x1 - real(ix)
      if(newcoef) then
         do k=4,1,-1
            cc0(k) = 0.
            DO I=IBGN, IEND
               J = I - IX + 2
               cc0(k) = cc0(k) + cf(i)*c0(k,j)
            enddo
         enddo

         do k=1,3
            cc1(k) = 0.
            DO I=IBGN, IEND
               J = I - IX + 2
               cc1(k) = cc1(k) + cf(i)*c1(k,j)
            enddo
         enddo
         cf1 = cf(1)
         ix0 = ix
      end if

      v  = cc0(4) + s*(cc0(3)+s*(cc0(2)+s*cc0(1)))
      db = (cc1(3)+s*(cc1(2)+s*cc1(1)))/dx

      IF (YLOG) THEN
         v   = EXP(B*LOG(10.D0))
         dv  = dv*EXP(B*LOG(10.D0))
      END IF
      IF (XLOG) THEN
         dv  = dv/XX
      END IF

      RETURN
      END
