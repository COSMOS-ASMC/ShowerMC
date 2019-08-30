        integer Ncomp,  Nev
        integer maxComp, idxmax
	parameter (maxComp = 300, idxmax=500)
        integer Ng(maxComp), Nec(maxComp)
        real*8  ElossT(maxComp), Eloss(idxmax, maxComp)
        common /Ztestc/ Eloss, ElossT, Ng, Nec, Ncomp, Nev

