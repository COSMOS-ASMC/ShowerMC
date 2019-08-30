      integer MaxIntMdls
      parameter(MaxIntMdls=8) 
!                    variables in Event generator
        type(ptcl):: Pjlab,     ! projectile in lab.
     *                Pjtatr,    ! projectile in target at rest system
     *                Pjcms,     ! projectile in cms
     *                Rpjlab,    ! recoiled proj. in lab
     *                Rpjtatr,   ! recoiled proj. in target at rest system
     *                Rpjcms     ! recoiled proj. in cms
!
        type(ptcl):: Tglab,     ! target in lab.
     *                Tgpatr,    ! target in projectile at rest system
     *                Tgcms,     ! target  in cms
     *                Rtglab,    ! recoiled target in lab
     *                Rtgpatr,   ! recoiled target in proj. at rest system
     *                Rtgcms     ! recoiled target in cms
!
        type(ptcl):: Cmsp,      ! cms particle formed by proj. and targe 
!                                  in lab.
     *                Missingp   ! missing mass particle (rest of
!                                  two recoiled leadings) in cms.
        real*8  Efrs,      ! effective roots 
     *          Powerexp,  ! pw2 in x**pw2 exp(-x)dx exp. type Pt distri.
     *          Powerp,    ! pw in x/(x +p0)**pw dx power type Pt distri.
     *          Ptnorm,    ! average Pt for normalization (pion's at this E)
     *          Probpower, ! Probability of power type pt.
     *          Avncharged ! Average number of charge particles 

        real*8  InteErg(MaxIntMdls) ! energy bound for interaction model.
        character*16 ModelList(MaxIntMdls) ! corresponding model name.
        real*8  InteErg2(MaxIntMdls) ! energy bound for XsecModel
        character*16 ModelList2(MaxIntMdls) ! and corresponding model name

!                             
        character*16 ActiveMdl  ! interaction model to be applied to the
        integer::ActiveMdlNo    ! ActiveMdel index
        character*16 ActiveMdl2  ! interaction model to be applied to the
                           !    current particle
        integer::ActiveMdlNo2   ! AtciveNdl2 index
        
        integer  Nnnb,    ! # of nn~
     *           Nkaon,   ! # of kaons
     *           Nddb,    ! # of DD~
     *           Npic,    ! # of pi+, pi-
     *           Nk0,     ! # of k0
     *           Nkch,    ! # of k+, k- ( Nkaon = Nkch + Nk0)
     *           Npi0,    ! # of pi 0
     *           Nch,     ! # of all charged secondaries
     *           Neta,    ! # of eta
     *           NoOfMdls, ! # of interaction models being used.
     *           NoOfMdls2 ! # of Xsec models being used.

        common /Zevhnv/  Pjlab, Pjtatr, Pjcms, Rpjlab, Rpjtatr,
     *                   Rpjcms, Tglab, Tgpatr, Tgcms, Rtglab,
     *                   Rtgpatr, Rtgcms, Cmsp, Missingp, 
     *                   InteErg, InteErg2,
     *                   Efrs, Powerexp, Powerp, Ptnorm, Probpower,
     *                   Avncharged, Nch,
     *                   Nnnb, Nkaon, Nddb, Npic, Nk0, Nkch, Npi0,
     *                   Neta, NoOfMdls, NoOfMdls2,
     *                   ActiveMdlNo, ActiveMdlNo2
       common /Zevhnvc/ ModelList, ActiveMdl, ModelList2, ActiveMdl2

       save /Zevhnv/, /Zevhnvc/  




