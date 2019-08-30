!            hadronic collision parameters 
!	(->  ---------------------------------------------

        character*64  IntModel !1 Interaction model description. Usage was changed from v6.0. 
                      !    One may list  code name and upper energy bound  for the code.\newline
                      ! E.g.  IntModel = '"dpmjet3"' ; to specify the dpmjet3 in the entire energy region
                      !       IntModel = '"dpmjet3" 100 "qgsjet2" to specify dpmjet3 at $<$ 100 GeV and qgsjetII 
                      !        at E$>$100 GeV \newline
                      !       IntModel = '"nucrin" 5 "fritiof1.6" 500 "adhoc" to specify  Nucrin,
                      !               fritiof1.6, and ad-hoc model in the respective energy region. This 
                      !       corresponds to the old IntModel='int1'. \newline
                      !       IntModel = '"nucrin" 5 "fritiof1.6" 10 "fritiof7.02"  and \newline
                      !       IntModel = '"dpmjet3"' \newline
                      !        are most promissing models that fit the observed 
                      !       data (muons and gamma rays) for which the primary is well known by
                      !       BESS and AMS observations ($<$ 100 GeV).
        character*64  XsecModel !1 Xsection model description.  Noarmally should not be given.
                      ! Defaul is to use the hadronic xsection given by each active model
                      ! fixed by IntModel. However, for some experimental purposes,
                      ! one may employ x-section by other interaction model.
                      !      e.g.
                      !     IntModel ='"phits" 2 "dpmjet" '
                      !     XsecModel='"Phits" 2 "dpmjet" 80 "qgsjet2"'
                      !  is one example. Default is blank and is replaced by IntModel

        character*100 InclusiveFile !2 The path to the inclusive data file, xdist.d. Default is
                       !   "../Contrib/Inclusive/xdist.d"
        real*8 SucPw   !2  In the 2nd, 3rd,.. collision of a nucleon inside a nucleus, the collision is 
                       !  made to be more elastic than normal one. The leading particle spectrum is
                       ! sampled from x**SucPw dx. SucPw should be in 1 to 2.
        real*8  Eta2Pi0 !2  eta/pi0 ratio. this is used to see the effect due to non-decay of pi0
                      ! at very high energies. Only source of h.e gamma can be eta and LPM may work
                      ! for them. default is 0.2
        integer MulLow !2 if 1, QCD predicted multiplicity law is used in the adhoc model else UA5
                       !  parametalization is used. Default is 1. (from v5), 
                       !  0.6135exp(23/18sqrt(2log(roots/0.3))) is QCD jet prediction. 
                       !  7.2roots**0.254 -7 is UA5 data. The branch point is set at roots = 900 GeV. 
                       !  (I have adjusted 0.6135 so that 900 GeV is the b.p)
        integer LundPara !2 To control Lund program. LundPara(1) is set to kfr(7); kfr(7)=1 is for Frititof
                       !  hard scattering. 2 is for Pythia H.S. 2 gives higher multiplicity but shape is 
                       !  strange.  Default is 1. LundPara(2) is set to kfr(12): 1 by for OPAL hard scattering
                       !   parameterization. 2 by DELPHI. Default is 2. (2 gives bit higher PT). LundPara(3)
                       !   $>$ 0 $\Rightarrow$    Pythia message will appear. LundPara(4) $>$ 0 Fritiof
                       !  message; both on ErrorOut. LundPara(5) =0 $\Rightarrow$    All kaons collisions
                       !  are treated as pi- in Fritiof, else they are treated by adhoc model as they are.
!                 below: Obsolete
        integer TotXSopt !2 option for total x-section.  1.  PDG  3. TOTEM (COMPETE fitting)
                         !   2. between 1 and 3.  (1 is lowest 3 is highest X-sec.) 
                         !   Diff. becomes gradually seen at  roots > few x 10^2 GeV
                         !   This for p-p case. Default is 2.

                         !   For other collision types, 1 is used.
                         !   If the current interaction model supplies the inelastic 
                         !   cross-section,  it is used without referring to this. 
                         !    However, see  XsecModel.
        integer SxAbySxpOpt !2  For nucleus target with mass # A, cross section is converted
                        !    from the one for proton target.  This option fixes which
                       !   SxA/sxp table is used. This is used when the current interaction
                       !   model dose not supply the cross-section. (See XsecModel too).
                       !   1.  use table derived from  QGSJET-II-04  (default)
                       !   2.  use table derived from dpjmet3
                       !   3.  use table dervied from EPOS  (but can be used for A<=64).
                       !        (as of 2013/Jun).
                       !   However, for small cross-sections (such as gamma-A, nu-A, old
                       !   cxp2xA routine is used).

        real*8  Cepic0  !2  Obsolete 
        real*8  Cekaon  !2  Obsolete
        integer SucInt  !2  The number of successive interactions inside A is affected by this parameter.\newline
                        !   If $0\rightarrow$ old formula (before uv3.500) is used, which give rather         
                        !   smaller number ($<Nsuc>$ in Air = 1.7 for 30 mb pp), \newline
                        !   if $1\rightarrow$ new one $<Nsuc>=2.2$ for 30 mb pp). \newline
                        !   Default is 0 (from V5.00 again).
        real*8  Ceneuc  !2 \verb| p -> n ; n-> p; p~->n~; n~->p~| prob.
        real*8  Mudirp  !2 \verb| DD~| \# enhancement factor. D is only for prompt muon.
     	real*8  Kpilog  !2  K$_{ch}/\pi_{ch}$=(Kpilog*log(ss+.069)+Kpicns)*exp(-8/s') where ss(GeV**2)= 
                        !   effective s. s'(GeV**2)=s - 4.63.  See also Kpicns.
        real*8  Kpicns  !2  See Kpilog. 0.077
        real*8  Elund   !2   obsolete (from v6.0)
!                  old def
        real*8  Elund2  !2  obsolete (from v6.0)
!              old def
         real*8 Elund3   !2 obsolete  (from v6.0).
          real*8  Efermi  !2  If Kinetic E $<$ Efermi, Fermi Momentum is considered for 
                         ! Nucleus target.
         integer:: DoNPadjust !2 0(Default):  0--> original is used
                    ! 1--> do adustement of  # of n,p  for phits model
                    !  when projectile is  p or n and target A > 10.
         integer:: dpmRareDecay !2 1(Default):  Rare ptcls in dpmjet3 is
                     ! forced to decay. (see cdpmRareDecay in cdpmjet3.f)
                     ! 0--> not decayed
         integer:: K0sSemiLD !2  12(Default).  
                !  semileptonic decay channel of K0s is
                !  1) pi+e + nue   B.R 7x10^-4
                !  2) pi+mu+ numu  B.R 6.4x10^-4
                ! if 0, these are neglected.
                ! if 1, only 1) is taken into account
                ! if 2, only 2) is //
                ! if 12, both are taken into account.  
!	<-)     -----------------------------
!            next are not input parameters.
        integer nmdls
!        parameter (nmdls = 12)
        parameter (nmdls = 13)
        character*16 RegMdls(nmdls)
        real*8 smallxs, largexs
        parameter (smallxs=1.d-10, largexs=1.d10)

        common /Zevhnp/ Cepic0, Cekaon, Ceneuc, Mudirp, Kpilog, Eta2Pi0,
     *         Kpicns,  Efermi, Elund, Elund2, Elund3, SucPw, MulLow, 
     *         LundPara(10), SucInt, DoNPadjust, dpmRareDecay,
     *         TotXSopt, SxAbySxpOpt, K0sSemiLD
!!     *         LundPara(10), SucInt
        common /Zevhnc/RegMdls, InclusiveFile, IntModel, XsecModel


