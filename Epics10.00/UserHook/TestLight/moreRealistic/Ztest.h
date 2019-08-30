        integer Ncomp,  Nev
        integer maxComp, idxmax
        parameter (maxComp = 20000)
        integer Ng(maxComp), Nec(maxComp)

        real*8  ElossT(maxComp)
        common /Ztestc/  ElossT, Ng, Nec, Ncomp, Nev


