        logical scorec,  lcorec, Molier, usees
	real*8 escb, escp, x0, emass, ecrit, tprob, Es, x0ing
	real*8 cconst, x0sq, constm, zchrg, amassn, ghmfp, ashad
        real*8 zbya,  emass2

        integer IncGp

!                          ***************
                           common /Zelmag/
!                          ***************
     1   escb,  escp,  x0,  emass,  ecrit,  tprob,
     2   Es,  x0ing, cconst, x0sq, constm, zchrg,
     3   amassn,  ghmfp, ashad, zbya, emass2,
     4   scorec,  lcorec,  Molier, usees,
     5   IncGp

!

!
!
!
! scorec:  if t, screening correction is applied
! lcorec:  if t, landau effect is taken into account
!   escb:  screening correction for brems is used below this energy
!   escp:  screening correction for pair cre is used below this energy
!     x0:  radiation length of lead in cm
!  emass:  electron mass in gev
!  ecrit:  criticl energy in lead.  moller and bhaba scattergin
!          contribution is inclued in this energy loss.
!  tprob:  total probability.  set when call bremst, compt or pairt
!     es:  scattergin energy in gev.  21.e-3
! Molier:  if t, Molier theory is used for scattering angle sampling
!  x0ing:  x0 expressed in g/cm**2.   5.82
! cconst:  3/8*thoms*z/a*x0ing for compton scattering.  to convert
!          total compton cross-section in unit of thomson into
!          / (radiation length).
! constm:  to convert moller scat. cross-section to /r.l
!          .3 * z/a*x0ing.   .3 = 2pi*re**2*n
!          used also in positron annihilation
! zchrg:   z of media
! amassn:  mass no. of media
! zbya:    z/a
! ashad:   xs(ga)=xs(gp)*ashad;  ashad=amassn*.75
! IncGp:   if 0 no hardronic process
!             1 vector dominace process is included
! ghmfp:   ghmfp/xsection (in mb) ---> m.f.p in r.l
! usees:   if t, escb and escp are used as boundary energy below
!          which partical screening correction is used else
!          constant values in each sampling routine are used.
! emass2:  2*emass
! StoppingPw:  Specify the formula # to be used for low energy heavy
!              ion dE/dx
!              1--> Effective charge method with complex terms
!              2--> //                      with simple terms by
!                   Pierce and Blann; Almost no diff. betw. 1 and 2
!                   Default is 1
!    ---- end of Zelmag ----
!
