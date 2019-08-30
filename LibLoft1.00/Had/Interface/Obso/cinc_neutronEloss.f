!              move incident to a:
            np=1
            a(1) = pj
!              
            alpha =( (ia-1.0)/(ia+1.0) )**2
!              Elastic scattering: E1--> E2
!             E2/E1 = (A**2 + 2Acos + 1)/(A+1)**2 = 1/2 ((alpha+1) + (1-alpha)*cos)
!             we treat average loss. gzai=log(E1/E2) = 2/( A + 2/3)
!             i.e,  E2 = E1*exp(-gzai).  cos angle can be fixed by E2/E1
!             after that, cos is converted to lab angle by
!             cosL = (1+Acos)/sqrt(A^2+2Acos+1)
            gzai = 2./(ia + 0.6666)
            E2 = ek*exp(-gzai)
            a(1).fm.p(4) = a(1).mass + E2
            cosz =( E2/ek *2. -(alpha + 1) )/(1.0-alpha) ! cm angle
!                  E2/ek*2 < 2  so -1-alpha  <E2/ek*2 - alpha -1 <1-alpha
!             and  -(1+alpha)/(1-alpah) < cosz < 1
!             since 0<= alpha <1,     -1<cosz<1
!             in Lab;   cosL =(1+Acos)/sqrt( A^2+2Acos+1)
!                            = (1+Acos)/sqrt( (Acos+ 1)^2 +A^2(1-cos^2))<1
            cosz = (1.0 + ia*cosz)/sqrt(ia**2 + 2*ia*cosz+1.)  !lab angle

!              azimuthal angel
            sinz = sqrt(1.-cosz*cosz)
            call kcossn(cs, sn)
            p = sqrt( E2 * a(1).fm.p(4) )
            a(1).fm.p(1) = p*cs*sinz
            a(1).fm.p(2) = p*sn*sinz
            a(1).fm.p(3) = p*cosz
