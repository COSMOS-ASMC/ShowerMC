
*-- Author :    D. HECK IK FZK KARLSRUHE       10/09/1998
C=======================================================================

      BLOCK DATA NEXDAT

C-----------------------------------------------------------------------

C  EPOS DATA INITIALIZATION
C
C  INITIALIZES DATA FOR EPOS LINK

C-----------------------------------------------------------------------


C  CONVERTS CORSIKA PARTICLE CODE TO EPOS/NEXUS PARTICLE CODE
      DATA IDTABL/
     *   10,  -12,   12,    0,  -14,   14,  110,  120, -120,  -20, !  10
     *  130, -130, 1220, 1120,-1120,   20,  220, 2130, 1130, 1230, !  20
     * 2230, 1330, 2330, 3331,-1220,-2130,-1130,-1230,-2230,-1330, !  30
     *-2330,-3331,    0,    0,    0,    0,    0,    0,    0,    0, !  40
     *    0,    0,    0,    0,    0,    0,    0,    0,    0,  221, !  50
     *  111,  121, -121, 1111, 1121, 1221, 2221,-1111,-1121,-1221, !  60
     *-2221,  231,  131, -131, -231,   11,  -11,   13,  -13,    0, !  70
     *  220,  220,  220,  220,    0,             25*0            ,
     *             15*0            , -140, -240,  240,  140,  340, ! 120
     * -340,  440, -141, -241,  241,  141,  341, -341,    0,  441, ! 130
     *  -16,   16,   15,  -15,    0,    0, 2140, 3140, 3240, 1140, ! 140
     * 1240, 2240, 1340, 2340, 3340,    0,    0,    0,-2140,-3140, ! 150
     *-3240,-1140,-1240,-2240,-1340,-2340,-3340,    0,    0,    0, ! 160
     * 1141, 1241, 2241,    0,    0,    0,    0,    0,    0,    0, ! 170
     *-1141,-1241,-2241,                      27*0               /

      END

*-- Author :    D. HECK IK FZK KARLSRUHE       10/09/1998
C=======================================================================

      SUBROUTINE NEXINI

C-----------------------------------------------------------------------

C  EPOS INITIALIZATION
C
C  FIRST INITIALIZATION OF EPOS ARRAYS AND PARAMETERS.

C  THIS SUBROUTINE IS CALLED FROM START.
C-----------------------------------------------------------------------

      IMPLICIT NONE
      COMMON /CRAIR/   COMPOS,PROBTA,AVERAW,AVOGDR
      DOUBLE PRECISION COMPOS(3),PROBTA(3),AVERAW,AVOGDR
      COMMON /CRDPMFLG/NFLAIN,NFLDIF,NFLPI0,NFLCHE,NFLPIF,NFRAGM
      INTEGER          NFLAIN,NFLDIF,NFLPI0,NFLCHE,NFLPIF,NFRAGM
      COMMON /CRPAM/   PAMA,SIGNUM,RESTMS,DECTIM
      DOUBLE PRECISION PAMA(6000),SIGNUM(6000),RESTMS(6000),
     *                 DECTIM(200)
      COMMON /CRRANDPA/RD,FAC,U1,U2,NSEQ,ISEED,KNOR
      DOUBLE PRECISION RD(3000),FAC,U1,U2
      INTEGER          ISEED(3,10),NSEQ
      LOGICAL          KNOR
      COMMON /CRRUNPAR/FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,
     *                 STEPFC,
# 4137 "corsika.h"
     *                 NRRUN,NSHOW,MPATAP,MONIIN,
     *                 MONIOU,MDEBUG,NUCNUC,MTABOUT,MLONGOUT,
     *                 ISEED1I,
# 4158 "corsika.h"
     *                 LSTCK,

     *                 LSTCK1,LSTCK2,

     *                 ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     *                 N1STTR,MDBASE,

     *                 DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     *                 FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     *                 FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN

     *                 ,FOUTFILE,IFINAM
# 4199 "corsika.h"
      COMMON /CRRUNPAC/DATDIR,DSN,DSNTAB,DSNLONG,HOST,USER
# 4215 "corsika.h"
     *                 ,FILOUT
# 4226 "corsika.h"
      DOUBLE PRECISION FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,STEPFC

      INTEGER          NRRUN,NSHOW,MPATAP,MONIIN,MONIOU,MDEBUG,NUCNUC,
     *                 ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     *                 N1STTR,MDBASE,MTABOUT,MLONGOUT,ISEED1I(3)
# 4257 "corsika.h"
      INTEGER          LSTCK

     *                ,LSTCK1,LSTCK2
# 4268 "corsika.h"
      CHARACTER*132    FILOUT

      CHARACTER*255    DSN,DSNTAB,DSNLONG
      CHARACTER*132    DATDIR
      CHARACTER*20     HOST,USER
# 4288 "corsika.h"
      LOGICAL          DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     *                 FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     *                 FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN
# 4301 "corsika.h"
      LOGICAL          FOUTFILE
      INTEGER          IFINAM
     

      COMMON /CRNEXPAR/NEXPRM,NNPARM,IRAND
      INTEGER          NEXPRM,NNPARM,IRAND(3)

      COMMON /CRNEXUS/ ISH0N,IVERNX,NEXVER,FNEXUS,FNEXSG
      INTEGER          ISH0N,IVERNX,NEXVER
      LOGICAL          FNEXUS,FNEXSG

C--------------------------- EPOS COMMON ------------------------------
      integer ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr
      real pnll,ptq,exmass,cutmss,wproj,wtarg
      common/hadr1/pnll,ptq,exmass,cutmss,wproj,wtarg
      integer idprojin,idtargin,irdmpr
     &             ,isoproj,isotarg
      real rexdifi,rexndii
      common/hadr25/idprojin,idtargin,rexdifi(4),rexndii(4),irdmpr
     &             ,isoproj,isotarg
      integer intpol,isigma,iomega,isetcs
      common/hadr6/intpol,isigma,iomega,isetcs
      integer iorsce,iorsdf,iorshh,ionudi
      common/cjinti/iorsce,iorsdf,iorshh,ionudi
      double precision seedi,seedj,seedj2,seedc
      integer iseqini,iseqsim
      common/cseed/seedi,seedj,seedj2,seedc,iseqini,iseqsim
      integer ishevt,ixtau,iwseed,jwseed,ixgeometry
      common/prnt3/ishevt,ixtau,iwseed,jwseed,ixgeometry
      integer istore,istmax,irescl,ntrymx,nclean,iopdg,ioidch
      real gaumx
      common/othe1/istore,istmax,gaumx,irescl,ntrymx,nclean,iopdg,ioidch
      integer      infragm
      common/nucl6/infragm
# 60540 "corsika.F"
      real egymin,egymax,elab,ecms,ekin
      common/enrgy/egymin,egymax,elab,ecms,ekin
      integer iversn,iverso
      common/versn/iversn,iverso
      integer iomodl,idproj,idtarg,laproj,maproj,latarg,matarg
      real wexcit,core,fctrmx
      common/hadr2/iomodl,idproj,idtarg,wexcit
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      integer icinpu
      real engy,elepti,elepto,angmue
      common/lept1/engy,elepti,elepto,angmue,icinpu
      real airznxs,airanxs,airwnxs,airavznxs,airavanxs
      common/nxsair/airznxs(3),airanxs(3),airwnxs(3)
     &             ,airavznxs,airavanxs

      INTEGER i
C-----------------------------------------------------------------------

      IF ( DEBUG ) WRITE(MDEBUG,*) 'NEXINI:'

      CALL aaset( 0 )
C  SET SOME PARAMETERS FOR COUPLING OF EPOS/NEXUS WITH CORSIKA
      ifmt = MONIOU
      ifch = MDEBUG
      ifop = NEXPRM
      CALL atitle

C  INITIALIZING FOR EPOS
      seedi  = ISEED(1,1) !dummy intialization
      seedj  = ISEED(1,1) !dummy intialization
      iwseed = 0
      egymin = 6.         ! MINIMUM ENERGY CORREPSONDS WITH ELAB= 18 GeV
      egymax = 2E6        ! MAXIMUM ENERGY CORRESPONDS WITH ELAB=2000EeV
      call LHCparameters  !LHC tune
      isigma = 0           ! do not show cross section
      isetcs=3    !option to obtain pomeron parameters
                  ! 0.....determine parameters but do not use Kfit
                  ! 1.....determine parameters and use Kfit
                ! else..get from table
                !         should be sufficiently detailed
                !          say iclegy1=1,iclegy2=99
                !         table is always done, more or less detailed!!!
                !and option to use cross section tables
                  ! 2....tabulation
                ! 3....simulation
      ionudi=1    !include quasi elastic events but strict calculation of xs
      iorsce=0      !color exchange turned on(1) or off(0)
      iorsdf=3      !droplet formation turned on(>0) or off(0)
      iorshh=0      !other hadron-hadron int. turned on(1) or off(0)

      istore = 0           ! DO NOT STORE EVENTS ON zzz.data FILE
CC the following commented statements are read by NEXPRM (see DATAC)
CC    ndecay = 1111110     ! let only resonances decay
CC    nrnody = nrnody + 1  ! number of nody to be set
CC    nody(nrnody) = 220   ! no decay for eta particles
CC    iappl  = 1           !
C  READ PARAMETERS SPECIFIED BY INPUT WITH KEY WORD 'EPOPAR' OR 'NEXPAR'
      CALL aread

      IF(NFRAGM.GE.3)THEN
        infragm=0
      ELSE
        infragm=NFRAGM
      ENDIF

      CLOSE( NEXPRM )
      airanxs(1)=14.
      airanxs(2)=16.
      airanxs(3)=40.
      airznxs(1)=7.
      airznxs(2)=8.
      airznxs(3)=18.
      airavznxs=0.
      do i=1,3
        airwnxs(i)=sngl(COMPOS(i))
        airavznxs=airavznxs+airwnxs(i)*airznxs(i)
      enddo
      airavanxs=sngl(AVERAW)
      NEXVER = iversn

C  DUMMY INITIALIZATIONS (FOR FILLING CROSS-SECTION TABLE)
      IF ( FNEXSG ) THEN

        idprojin = 1120
        idtargin = 1120

        maproj = 56
        laproj = 28
        matarg = 14
        latarg = 1
        pnll   = 200.
        engy   = -1.
        CALL ainit
      ENDIF

      IF ( DEBUG ) WRITE(MDEBUG,*) 'NEXINI: EXIT'
      WRITE(MDEBUG,*) ' '

      RETURN
      END

*-- Author :    D. HECK IK FZK KARLSRUHE       10/09/1998
C=======================================================================

      SUBROUTINE NEXLNK

C-----------------------------------------------------------------------
C  (EPOS/)NEX(US) L(I)NK (TO CORSIKA)
C
C  LINKS EPOS/NEXUS PACKAGE TO CORSIKA, NEEDS FIRST CALL OF NEXINI.
C  THIS SUBROUTINE IS CALLED FROM SDPM.
C-----------------------------------------------------------------------

      IMPLICIT NONE

      COMMON /CRPAM/   PAMA,SIGNUM,RESTMS,DECTIM
      DOUBLE PRECISION PAMA(6000),SIGNUM(6000),RESTMS(6000),
     *                 DECTIM(200)

      COMMON /CRPARPAR/CURPAR,SECPAR,PRMPAR,OUTPAR,C,

     *                 E00,E00PN,PTOT0,PTOT0N,THICKH,ITYPE,LEVL

# 3932 "corsika.h"
      DOUBLE PRECISION CURPAR(0:16),SECPAR(0:16),PRMPAR(0:16),
     *                 OUTPAR(0:16),

     *                 C(50),E00,E00PN,PTOT0,PTOT0N,THICKH
      INTEGER          ITYPE,LEVL

      DOUBLE PRECISION GAMMA,COSTHE,PHIX,PHIY,H,T,X,Y,CHI,BETA,GCM,ECM
# 3954 "corsika.h"
      EQUIVALENCE      (CURPAR(1), GAMMA ), (CURPAR(2), COSTHE),
     *                 (CURPAR(3), PHIX  ), (CURPAR(4), PHIY  ),
     *                 (CURPAR(5), H     ), (CURPAR(6), T     ),
     *                 (CURPAR(7), X     ), (CURPAR(8), Y     ),
     *                 (CURPAR(9), CHI   ), (CURPAR(10),BETA  ),
     *                 (CURPAR(11),GCM   ), (CURPAR(12),ECM   )

      COMMON /CRRANDPA/RD,FAC,U1,U2,NSEQ,ISEED,KNOR
      DOUBLE PRECISION RD(3000),FAC,U1,U2
      INTEGER          ISEED(3,10),NSEQ
      LOGICAL          KNOR

      COMMON /CRREST/  CONTNE,TAR,LT
      DOUBLE PRECISION CONTNE(3),TAR
      INTEGER          LT

      COMMON /CRRUNPAR/FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,
     *                 STEPFC,
# 4137 "corsika.h"
     *                 NRRUN,NSHOW,MPATAP,MONIIN,
     *                 MONIOU,MDEBUG,NUCNUC,MTABOUT,MLONGOUT,
     *                 ISEED1I,
# 4158 "corsika.h"
     *                 LSTCK,

     *                 LSTCK1,LSTCK2,

     *                 ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     *                 N1STTR,MDBASE,

     *                 DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     *                 FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     *                 FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN

     *                 ,FOUTFILE,IFINAM
# 4199 "corsika.h"
      COMMON /CRRUNPAC/DATDIR,DSN,DSNTAB,DSNLONG,HOST,USER
# 4215 "corsika.h"
     *                 ,FILOUT
# 4226 "corsika.h"
      DOUBLE PRECISION FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,STEPFC

      INTEGER          NRRUN,NSHOW,MPATAP,MONIIN,MONIOU,MDEBUG,NUCNUC,
     *                 ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     *                 N1STTR,MDBASE,MTABOUT,MLONGOUT,ISEED1I(3)
# 4257 "corsika.h"
      INTEGER          LSTCK

     *                ,LSTCK1,LSTCK2
# 4268 "corsika.h"
      CHARACTER*132    FILOUT

      CHARACTER*255    DSN,DSNTAB,DSNLONG
      CHARACTER*132    DATDIR
      CHARACTER*20     HOST,USER
# 4288 "corsika.h"
      LOGICAL          DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     *                 FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     *                 FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN
# 4301 "corsika.h"
      LOGICAL          FOUTFILE
      INTEGER          IFINAM
    

      COMMON /CRNEXLIN/ IDTABL
      INTEGER           IDTABL(200)

      COMMON /CRNEXPAR/NEXPRM,NNPARM,IRAND
      INTEGER          NEXPRM,NNPARM,IRAND(3)

      COMMON /CRNEXUS/ ISH0N,IVERNX,NEXVER,FNEXUS,FNEXSG
      INTEGER          ISH0N,IVERNX,NEXVER
      LOGICAL          FNEXUS,FNEXSG


C--------------------------- EPOS COMMON ------------------------------
      real pnll,ptq,exmass,cutmss,wproj,wtarg
      common/hadr1/pnll,ptq,exmass,cutmss,wproj,wtarg
      integer nevt,kolevt,npjevt,minfra,maxfra,nglevt,koievt,kohevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
      real phievt,bimevt,pmxevt,egyevt,xbjevt,qsqevt,zppevt,zptevt
      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra,kohevt
      integer idprojin,idtargin,irdmpr
     &             ,isoproj,isotarg
      real rexdifi,rexndii
      common/hadr25/idprojin,idtargin,rexdifi(4),rexndii(4),irdmpr
     &             ,isoproj,isotarg
      double precision seedi,seedj,seedj2,seedc
      integer iseqini,iseqsim
      common/cseed/seedi,seedj,seedj2,seedc,iseqini,iseqsim
# 60716 "corsika.F"
      integer iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      integer iomodl,idproj,idtarg,laproj,maproj,latarg,matarg
      real wexcit,core,fctrmx
      common/hadr2/iomodl,idproj,idtarg,wexcit
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      integer icinpu
      real engy,elepti,elepto,angmue
      common/lept1/engy,elepti,elepto,angmue,icinpu
      real egymin,egymax,elab,ecms,ekin
      common/enrgy/egymin,egymax,elab,ecms,ekin
C-----------------------------------------------------------------------
      DOUBLE PRECISION ELABN
      INTEGER J,NNEUT,NPROT
      SAVE
C-----------------------------------------------------------------------

      IF ( DEBUG ) THEN
        WRITE(MDEBUG,*) 'NEXLNK: TAR',SNGL(TAR)
        ish  = ISH0N
C  RANDOM GENERATOR STATUS (SEQUENCE L=1) AT BEGINNING OF EVENT
        CALL RMMAQD( ISEED(1,1),1,'R' )
        IRAND(1) = ISEED(1,1)
C  NUMBER OF CALLS
        IRAND(2) = ISEED(2,1)
C  NUMBER OF BILLIONS
        IRAND(3) = ISEED(3,1)
        WRITE(MDEBUG,158) (IRAND(J),J=1,3)
 158    FORMAT(' NEXLNK: RANDOM NUMBER GENERATOR AT BEGIN:'
     *        ,' SEQUENCE= 1 SEED= ',I9,' CALLS=',I9,
     *         ' BILLIONS=',I9)
ctp     ELSE
ctp       ish  = 0
      ENDIF
      seedc=ISEED(2,1)+1.D9*ISEED(3,1)
C  CALCULATE ENERGY IN LAB SYSTEM FOR ELASTICITY FOR VARIOUS PROJECTILES
      IF     ( ITYPE .EQ. 1 ) THEN
C  TREAT GAMMA PROJECTILES (FROM EGS)
        CALL RMMARD( RD,1,1 )
        IF ( RD(1) .LE. 0.5D0 ) THEN
          ITYPE = 7
        ELSE
          ITYPE = 17
        ENDIF
        ELABN = CURPAR(1)
        CURPAR(1) = ELABN / PAMA(ITYPE)
      ELSEIF ( ITYPE .LT. 200 ) THEN
C  TREAT ORDINARY PROJECTILES
        ELABN  = CURPAR(1) * PAMA(ITYPE)
      ELSE
C  TREAT NUCLEI PROJECTILES
        NPROT = MOD(ITYPE,100)
        NNEUT = ITYPE/100 - NPROT
        ELABN = CURPAR(1) * ( PAMA(14)*NPROT + PAMA(13)*NNEUT )
      ENDIF
C  SET TARGET PARAMETERS
      matarg = NINT( TAR )

      idtargin = 1120

      IF     ( TAR .EQ. 14.D0 ) THEN
        latarg = 7
      ELSEIF ( TAR .EQ. 16.D0 ) THEN
        latarg = 8
      ELSEIF ( TAR .EQ. 40.D0 ) THEN
        latarg = 18
# 60822 "corsika.F"
      ELSE
        WRITE(MONIOU,*) 'NEXLNK: UNDEFINED TARGET TAR=',SNGL(TAR)

      ENDIF

C  SET PROJECTILE PARAMETERS
      IF ( ITYPE .LT. 200 ) THEN

        idprojin = IDTABL(ITYPE)
        IF     ( idprojin .EQ. 20  .OR.  idprojin .EQ. -20 ) THEN
C  TREAT NEUTRAL KAONS  (K(0)S AND K(0)L)
          CALL RMMARD( RD,1,1 )
          IF ( RD(1) .LE. 0.5D0 ) THEN
            idprojin = 230
          ELSE
            idprojin = -230
          ENDIF
        ELSEIF ( idprojin .EQ. 2130 ) THEN
C  EPOS CANNOT TREAT LAMBDA, TAKE INSTEAD SIGMA(0))
          idprojin = 1230
        ELSEIF ( idprojin .EQ. -2130 ) THEN
C  EPOS CANNOT TREAT ANTI-LAMBDA, TAKE INSTEAD ANTI-SIGMA(0))
          idprojin = -1230
# 60868 "corsika.F"
        ENDIF
C  ALL OTHER PARTICLE CODES UNCHANGED
        laproj = -1
        maproj = 1
CDH2003 PNLL   = CURPAR(1)*PAMA(ITYPE)
        pnll   = PAMA(ITYPE) * BETA * GAMMA
      ELSE
C  PROJECTILE IS NUCLEUS

        idprojin = 1120

        laproj = MOD(ITYPE,100)
        maproj = ITYPE/100
CDH2003 PNLL   = CURPAR(1)*(PAMA(14)+PAMA(13))*0.5
        pnll   = 0.5 * (PAMA(14)+PAMA(13)) * BETA * GAMMA
      ENDIF
C  SET ENGY NEGATIVE TO FORCE CALCULATION IN LAB FRAME
      engy = -1.
      ecms = -1.
      elab = -1.
      ekin = -1.
C  INTIALIZE ENERGY AND PARTICLE DEPENDENT PORTION OF EPOS/NEXUS
C  AT THE FIRST CALL: READ ALSO DATA SETS
      CALL AINIT

      IF ( DEBUG ) WRITE(MDEBUG,*) 'NEXLNK: AEPOS IS NOW CALLED'
C  CALL TO EPOS
      CALL AEPOS( 1 )
      IF ( DEBUG ) WRITE(MDEBUG,*) 'NEXLNK: RETURN FROM AEPOS'

      IF ( NEVT .EQ. 1 ) CALL AFINAL
C  NOW BRING PARTICLES TO CORSIKA STACK
      CALL NSTORE( ELABN )
      IF ( DEBUG ) WRITE(MDEBUG,*) 'NEXLNK: (EXIT)'

      RETURN
      END

*-- Author :    D. HECK IK FZK KARLSRUHE       10/09/1998
C=======================================================================

      SUBROUTINE NEXSIG( ELAB,ICZ )

C-----------------------------------------------------------------------

C  EPOS SIG(MAS)
C
C  CALCULATES INELASTIC HADRON-AIR CROSS-SECTIONS FOR EPOS MODEL.

C  NUCLEUS-AIR CROSS-SECTIONS ARE DETERMINED BY P-P CROSS-SECTIONS AND
C  THE CORSIKA GLAUBER TABLES (SEE BOX2).
C  THIS SUBROUTINE IS CALLED FROM BOX2.
C  ARGUMENTS:
C   ELAB   =  LABORATORY ENERGY (GEV)
C   ICZ    =  HADRON TYPE: 1 = PION, 2 = NUCLEON, 3 = KAON
C-----------------------------------------------------------------------

      IMPLICIT NONE

      COMMON /CRPARPAR/CURPAR,SECPAR,PRMPAR,OUTPAR,C,

     *                 E00,E00PN,PTOT0,PTOT0N,THICKH,ITYPE,LEVL

# 3932 "corsika.h"
      DOUBLE PRECISION CURPAR(0:16),SECPAR(0:16),PRMPAR(0:16),
     *                 OUTPAR(0:16),

     *                 C(50),E00,E00PN,PTOT0,PTOT0N,THICKH
      INTEGER          ITYPE,LEVL


      COMMON /CRRUNPAR/FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,
     *                 STEPFC,
# 4137 "corsika.h"
     *                 NRRUN,NSHOW,MPATAP,MONIIN,
     *                 MONIOU,MDEBUG,NUCNUC,MTABOUT,MLONGOUT,
     *                 ISEED1I,
# 4158 "corsika.h"
     *                 LSTCK,

     *                 LSTCK1,LSTCK2,

     *                 ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     *                 N1STTR,MDBASE,

     *                 DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     *                 FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     *                 FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN

     *                 ,FOUTFILE,IFINAM
# 4199 "corsika.h"
      COMMON /CRRUNPAC/DATDIR,DSN,DSNTAB,DSNLONG,HOST,USER
# 4215 "corsika.h"
     *                 ,FILOUT
# 4226 "corsika.h"
      DOUBLE PRECISION FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,STEPFC

      INTEGER          NRRUN,NSHOW,MPATAP,MONIIN,MONIOU,MDEBUG,NUCNUC,
     *                 ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     *                 N1STTR,MDBASE,MTABOUT,MLONGOUT,ISEED1I(3)
# 4257 "corsika.h"
      INTEGER          LSTCK

     *                ,LSTCK1,LSTCK2
# 4268 "corsika.h"
      CHARACTER*132    FILOUT

      CHARACTER*255    DSN,DSNTAB,DSNLONG
      CHARACTER*132    DATDIR
      CHARACTER*20     HOST,USER
# 4288 "corsika.h"
      LOGICAL          DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     *                 FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     *                 FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN
# 4301 "corsika.h"
      LOGICAL          FOUTFILE
      INTEGER          IFINAM
# 4323 "corsika.h"

# 4333 "corsika.h"

      COMMON /CRSIGM/  SIGMA,SIGANN,SIGAIR,FRACTN,FRCTNO,SIGAIRS
      DOUBLE PRECISION SIGMA,SIGANN,SIGAIR,FRACTN,FRCTNO,SIGAIRS

    

      COMMON /CRNEXSGM/XFRACN,XFRANO,SIGNAIR,SIGNHN
      DOUBLE PRECISION XFRACN(12,3),XFRANO(12,3),SIGNAIR(12,3),
     *                 SIGNHN(12,3)

     


      DOUBLE PRECISION DELTAE,ELAB,SECT,WK(3),YE
      INTEGER          I,ICZ,JE
      SAVE
C-----------------------------------------------------------------------

      IF ( DEBUG ) WRITE(MDEBUG,*) 'NEXSIG: ELAB=',SNGL(ELAB),
     *                                 ' ICZ=',ICZ

C  DETERMINE ENERGY INTERVAL FOR INTERPOLATION
      YE = LOG10(ELAB)
      IF ( YE .LT. 1.D0 ) YE = 1.D0
      JE = INT( YE )
      IF ( JE .GT. 9 ) JE = 9
      DELTAE = YE - JE
      WK(3)  = DELTAE * (DELTAE-1.D0) * .5D0
      WK(1)  = 1.D0 - DELTAE  + WK(3)
      WK(2)  = DELTAE - 2.D0 * WK(3)

      IF     ( ICZ .LE.   3 ) THEN
C  FOR HADRON PROJECTILES
        SECT = 0.D0
        DO  I = 1, 3
          SECT = SECT + SIGNAIR(JE+I-1,ICZ)*WK(I)
        ENDDO
        SIGAIR = EXP( SECT )
        SECT   = 0.D0
        DO  I = 1, 3
          SECT = SECT + XFRACN(JE+I-1,ICZ)*WK(I)
        ENDDO
        FRACTN = EXP( SECT )
        SECT   = 0.D0
        DO  I = 1, 3
          SECT = SECT + XFRANO(JE+I-1,ICZ)*WK(I)
        ENDDO
        FRCTNO = EXP( SECT )
        SIGMA  = 0.D0

      ELSEIF ( ICZ .GE. 200 ) THEN
C  FOR NUCLEUS PROJECTILES DETERMINE ONLY NN CROSS-SECTION
        SIGAIR = 0.D0
        FRACTN = 0.D0
        FRCTNO = 0.D0
        SIGMA = 0.D0
        DO  I = 1, 3
          SIGMA = SIGMA + SIGNHN(JE+I-1,2)*WK(I)
        ENDDO
        SIGMA = EXP( SIGMA )

      ELSE

        WRITE(MONIOU,444) (CURPAR(I),I=0,9)
 444    FORMAT(' NEXSIG: CURPAR=',1P,10E11.3)

        WRITE(MONIOU,*) 'NEXSIG: ILLEGAL PROJECTILE TYP =',ICZ
        STOP
      ENDIF
      IF ( DEBUG ) WRITE(MDEBUG,*) 'NEXSIG: SIGMA=',SNGL(SIGMA),
     *                                  ' SIGAIR=',SNGL(SIGAIR)

      RETURN
      END

*-- Author :    D. HECK IK FZK KARLSRUHE       10/09/1998
C=======================================================================

      SUBROUTINE NEXSIGINI

C-----------------------------------------------------------------------
C  NEX(US) SIG(MAS) INI(TIALIZATION)
C
C  INITIALIZES INELASTIC CROSS-SECTION.
C  INTEGER 'ICZ' IS HADRON TYPE: 1 = PION, 2 = NUCLEON, 3 = KAON
C  THIS SUBROUTINE IS CALLED FROM START.
C-----------------------------------------------------------------------

      IMPLICIT NONE
      COMMON /CRAIR/   COMPOS,PROBTA,AVERAW,AVOGDR
      DOUBLE PRECISION COMPOS(3),PROBTA(3),AVERAW,AVOGDR

      COMMON /CRPAM/   PAMA,SIGNUM,RESTMS,DECTIM
      DOUBLE PRECISION PAMA(6000),SIGNUM(6000),RESTMS(6000),
     *                 DECTIM(200)

      COMMON /CRRUNPAR/FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,
     *                 STEPFC,
# 4137 "corsika.h"
     *                 NRRUN,NSHOW,MPATAP,MONIIN,
     *                 MONIOU,MDEBUG,NUCNUC,MTABOUT,MLONGOUT,
     *                 ISEED1I,
# 4158 "corsika.h"
     *                 LSTCK,

     *                 LSTCK1,LSTCK2,

     *                 ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     *                 N1STTR,MDBASE,

     *                 DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     *                 FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     *                 FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN

     *                 ,FOUTFILE,IFINAM
# 4199 "corsika.h"
      COMMON /CRRUNPAC/DATDIR,DSN,DSNTAB,DSNLONG,HOST,USER
# 4215 "corsika.h"
     *                 ,FILOUT
# 4226 "corsika.h"
      DOUBLE PRECISION FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,STEPFC

      INTEGER          NRRUN,NSHOW,MPATAP,MONIIN,MONIOU,MDEBUG,NUCNUC,
     *                 ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     *                 N1STTR,MDBASE,MTABOUT,MLONGOUT,ISEED1I(3)
# 4257 "corsika.h"
      INTEGER          LSTCK

     *                ,LSTCK1,LSTCK2
# 4268 "corsika.h"
      CHARACTER*132    FILOUT

      CHARACTER*255    DSN,DSNTAB,DSNLONG
      CHARACTER*132    DATDIR
      CHARACTER*20     HOST,USER
# 4288 "corsika.h"
      LOGICAL          DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     *                 FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     *                 FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN
# 4301 "corsika.h"
      LOGICAL          FOUTFILE
      INTEGER          IFINAM


      COMMON /CRNEXSGM/XFRACN,XFRANO,SIGNAIR,SIGNHN
      DOUBLE PRECISION XFRACN(12,3),XFRANO(12,3),SIGNAIR(12,3),
     *                 SIGNHN(12,3)


      COMMON /PSAR33/  ASECT,ASECTN
      REAL             ASECT(7,4,7)    !  neXus3
      REAL             ASECTN(7,7,7)   !  neXus3 for nucleus-nucleus

      DOUBLE PRECISION ATAR,ELAB,ENGY,EON,SECT(3),SECTN,SECTA,SECTO
      DOUBLE PRECISION WA(3),WK(3),YA,YE
      INTEGER          IAT,ICZ,IPROJ,JA,JE,JEE,K,M
      SAVE
C-----------------------------------------------------------------------

      IF ( DEBUG ) WRITE(MDEBUG,*) 'NEXSIGINI - START'

      IF ( DEBUG ) WRITE(MDEBUG,90)
 90   FORMAT('  ENERGY(LAB)    SIGAIR(PI)    SIGAIR(N)     SIGAIR(K)',
     *                '     SIGMA(pi-p)   SIGMA(p-p)    SIGMA(K-p) ')
      DO  100  JEE = 1, 12
C  LOOP 100 RUNS OVER ALL ENERGY VALUES
C  CALCULATE CM ENERGY 'ENGY'  FROM LAB ENERGY 'ELAB'
        ELAB = 10.D0**JEE

        DO  99  ICZ = 1, 3
C  HADRON PROJECTILES ICZ ARE: PI - 1, N - 2, K - 3
C  FIRST FOR NITROGEN TARGET
          IF     ( ICZ .EQ. 1 ) THEN
            IPROJ = 8
          ELSEIF ( ICZ .EQ. 2 ) THEN
            IPROJ = 14
          ELSEIF ( ICZ .EQ. 3 ) THEN
            IPROJ = 11
          ENDIF
          ATAR = PAMA(14)
          ENGY = SQRT( 2.D0*ELAB*ATAR +ATAR**2+PAMA(IPROJ)**2 )
          YE   = LOG10(MAX( 1.D0, ENGY/1.5D0 ))+1.D0
CC        IF ( DEBUG ) WRITE(MDEBUG,*) 'ENGY=',ENGY,' ELAB=',ELAB
          IF ( YE .LT. 1.D0 ) YE = 1.D0
          JE = MIN( 5, INT( YE ) )
          WK(2) = YE - JE
          WK(3) = WK(2) * (WK(2)-1.D0) * .5D0
          WK(1) = 1.D0 - WK(2) + WK(3)
          WK(2) = WK(2) - 2.D0*WK(3)

          IAT = 14
          YA  = IAT
          YA  = DLOG( YA ) / 0.69315D0 + 1.D0
          JA  = MIN( INT( YA ), 4 )
          WA(2) = YA - JA
          WA(3) = WA(2) * (WA(2)-1.D0) * .5D0
          WA(1) = 1.D0 - WA(2) + WA(3)
          WA(2) = WA(2) - 2.D0*WA(3)

          SECTN = 0.D0
          DO  K = 1, 3
            DO  M = 1, 3
              SECTN = SECTN + ASECT(JE+K-1,ICZ,JA+M-1)*WA(M)*WK(K)
            ENDDO
          ENDDO
CC        IF ( DEBUG ) WRITE(MDEBUG,*)
CC   *             'ICZ,JEE,JE=',ICZ,JEE,JE,' SECTN=',SECTN
          SECTN = EXP( SECTN )
C  THEN  FOR OXYGEN TARGET
          ATAR = PAMA(14)
          ENGY = SQRT( 2.D0*ELAB*ATAR+ATAR**2+PAMA(IPROJ)**2 )
          YE   = LOG10(MAX( 1.D0, ENGY/1.5D0 ))+1.D0
          IF ( YE .LT. 1.D0 ) YE = 1.D0
          JE = MIN( 5, INT( YE ) )
          WK(2) = YE - JE
          WK(3) = WK(2) * (WK(2)-1.D0) * .5D0
          WK(1) = 1.D0 - WK(2) + WK(3)
          WK(2) = WK(2) - 2.D0*WK(3)

          IAT = 16
          YA  = IAT
          YA  = DLOG( YA ) / 0.69315D0 + 1.D0
          JA  = MIN( INT( YA ), 4 )
          WA(2) = YA - JA
          WA(3) = WA(2) * (WA(2)-1.D0) * .5D0
          WA(1) = 1.D0 - WA(2) + WA(3)
          WA(2) = WA(2) - 2.D0*WA(3)

          SECTO = 0.D0
          DO  K = 1, 3
            DO  M = 1, 3
              SECTO = SECTO + ASECT(JE+K-1,ICZ,JA+M-1)*WA(M)*WK(K)
            ENDDO
          ENDDO
CC        IF ( DEBUG ) WRITE(MDEBUG,*)
CC   *             'ICZ,JEE,JE=',ICZ,JEE,JE,' SECTO=',SECTO
          SECTO = EXP( SECTO )
C  THEN  FOR ARGON TARGET
          ATAR = PAMA(14)
          ENGY = SQRT( 2.D0*ELAB*ATAR+ATAR**2+PAMA(IPROJ)**2 )
          YE   = LOG10(MAX( 1.D0, ENGY/1.5D0 ))+1.D0
          IF ( YE .LT. 1.D0 ) YE = 1.D0
          JE = MIN( 5, INT( YE ) )
          WK(2) = YE - JE
          WK(3) = WK(2) * (WK(2)-1.D0) * .5D0
          WK(1) = 1.D0 - WK(2) + WK(3)
          WK(2) = WK(2) - 2.D0*WK(3)

          IAT = 40
          YA  = IAT
          YA  = DLOG( YA ) / 0.69315D0 + 1.D0
          JA  = MIN( INT( YA ), 4 )
          WA(2) = YA - JA
          WA(3) = WA(2) * (WA(2)-1.D0) * .5D0
          WA(1) = 1.D0 - WA(2) + WA(3)
          WA(2) = WA(2) - 2.D0*WA(3)

          SECTA = 0.D0
          DO  K = 1, 3
            DO  M = 1, 3
              SECTA = SECTA+ASECT(JE+K-1,ICZ,JA+M-1)*WA(M)*WK(K)
            ENDDO
          ENDDO
CC        IF ( DEBUG ) WRITE(MDEBUG,*)
CC   *             'ICZ,JEE,JE=',ICZ,JEE,JE,' SECTA=',SECTA
          SECTA  = EXP( SECTA )
C  NOW TAKE THE COMPOSITION OF AIR TO CALCULATE AIR CROSS-SECTION
          SECT(ICZ)        = COMPOS(1)*SECTN
          XFRACN(JEE,ICZ)  = LOG( SECT(ICZ) )
          SECT(ICZ)        = SECT(ICZ) + COMPOS(2)*SECTO
          XFRANO(JEE,ICZ)  = LOG( SECT(ICZ) )
          SECT(ICZ)        = SECT(ICZ) + COMPOS(3)*SECTA
          SIGNAIR(JEE,ICZ) = LOG( SECT(ICZ) )
C  PION NUCLEON,  NUCLEON NUCLEON,  KAON NUCLEON CROSS-SECTION
          ENGY = SQRT( 2.D0*ELAB*PAMA(14)+PAMA(14)**2+PAMA(IPROJ)**2 )
          YE   = LOG10(MAX( 1.D0, ENGY/1.5D0 ))+1.D0
          IF ( YE .LT. 1.D0 ) YE = 1.D0
          JE = MIN( 5, INT( YE ) )
          WK(2) = YE - JE
          WK(3) = WK(2) * (WK(2)-1.D0) * .5D0
          WK(1) = 1.D0 - WK(2) + WK(3)
          WK(2) = WK(2) - 2.D0*WK(3)
          SIGNHN(JEE,ICZ) = 0.D0
          DO  K = 1, 3
            SIGNHN(JEE,ICZ) = SIGNHN(JEE,ICZ)+ASECT(JE+K-1,ICZ,1)*WK(K)
          ENDDO

  99    CONTINUE
        IF ( DEBUG ) THEN
          EON  = 10.D0**JEE
          WRITE(MDEBUG,103) EON,SECT(1),SECT(2),SECT(3),
     *        EXP(SIGNHN(JEE,1)),EXP(SIGNHN(JEE,2)),EXP(SIGNHN(JEE,3))
 103      FORMAT(' ',7G14.5)
        ENDIF
 100  CONTINUE
      IF ( DEBUG ) THEN
        WRITE(MDEBUG,*)
        WRITE(MDEBUG,*) 'NOW LOGARITHMS OF THE CROSS-SECTIONS'
        WRITE(MDEBUG,90)
        DO  JEE = 1, 12
          EON = 10.D0**JEE
          WRITE(MDEBUG,103) EON,SIGNAIR(JEE,1),SIGNAIR(JEE,2),
     *      SIGNAIR(JEE,3),SIGNHN(JEE,1),SIGNHN(JEE,2),SIGNHN(JEE,3)
        ENDDO
        WRITE(MDEBUG,*) 'NEXSIGINI - END'
      ENDIF

      RETURN
      END

*-- Author :    D. HECK IK FZK KARLSRUHE       10/09/1998
C=======================================================================

      SUBROUTINE NSTORE( ELABN )

C-----------------------------------------------------------------------

C  EPOS PARTICLES STORE (INTO CORSIKA STACK)

C
C  STORES EPOS/NEXUS OUTPUT PARTICLES INTO CORSIKA STACK.
C  THIS SUBROUTINE IS CALLED FROM NEXLNK.
C  ARGUMENT:
C   ELABN  = ENERGY/NUCLEON OF PROJECTILE (GEV)
C-----------------------------------------------------------------------

      IMPLICIT NONE

      COMMON /CRDPMFLG/NFLAIN,NFLDIF,NFLPI0,NFLCHE,NFLPIF,NFRAGM
      INTEGER          NFLAIN,NFLDIF,NFLPI0,NFLCHE,NFLPIF,NFRAGM

      COMMON /CRELADPM/ELMEAN,ELMEAA,IELDPM,IELDPA
      DOUBLE PRECISION ELMEAN(40),ELMEAA(40)
      INTEGER          IELDPM(40,13),IELDPA(40,13)

      COMMON /CRELASTY/ELAST
      DOUBLE PRECISION ELAST

      COMMON /CRISTA/  IFINET,IFINNU,IFINKA,IFINPI,IFINHY,IFINCM,IFINOT
      INTEGER          IFINET,IFINNU,IFINKA,IFINPI,IFINHY,IFINCM,IFINOT

      INTEGER          LNGMAX
      PARAMETER        (LNGMAX = 15000)
      COMMON /CRLONGI/ ADLONG,AELONG,APLONG,DLONG,ELONG,HLONG,PLONG,
     *                 SDLONG,SELONG,SPLONG,THSTEP,THSTPI,
     *                 LHEIGH,NSTEP,
     *                 LLONGI,FLGFIT
# 3643 "corsika.h"
      DOUBLE PRECISION ADLONG(0:LNGMAX,19),AELONG(0:LNGMAX,10),
     *                 APLONG(0:LNGMAX,10),DLONG(0:LNGMAX,19),
     *                 ELONG(0:LNGMAX,10),
     *                 HLONG(0:LNGMAX),PLONG(0:LNGMAX,10),
     *                 SDLONG(0:LNGMAX,19),SELONG(0:LNGMAX,10),
     *                 SPLONG(0:LNGMAX,10),THSTEP,THSTPI

      INTEGER          LHEIGH,NSTEP
      LOGICAL          LLONGI,FLGFIT
      COMMON /CRMULT/  EKINL,MSMM,MULTMA,MULTOT
      DOUBLE PRECISION EKINL
      INTEGER          MSMM,MULTMA(40,13),MULTOT(40,13)

      COMMON /CRPAM/   PAMA,SIGNUM,RESTMS,DECTIM
      DOUBLE PRECISION PAMA(6000),SIGNUM(6000),RESTMS(6000),
     *                 DECTIM(200)

      COMMON /CRPARPAR/CURPAR,SECPAR,PRMPAR,OUTPAR,C,

     *                 E00,E00PN,PTOT0,PTOT0N,THICKH,ITYPE,LEVL

# 3932 "corsika.h"
      DOUBLE PRECISION CURPAR(0:16),SECPAR(0:16),PRMPAR(0:16),
     *                 OUTPAR(0:16),

     *                 C(50),E00,E00PN,PTOT0,PTOT0N,THICKH
      INTEGER          ITYPE,LEVL

      DOUBLE PRECISION GAMMA,COSTHE,PHIX,PHIY,H,T,X,Y,CHI,BETA,GCM,ECM
# 3954 "corsika.h"
      EQUIVALENCE      (CURPAR(1), GAMMA ), (CURPAR(2), COSTHE),
     *                 (CURPAR(3), PHIX  ), (CURPAR(4), PHIY  ),
     *                 (CURPAR(5), H     ), (CURPAR(6), T     ),
     *                 (CURPAR(7), X     ), (CURPAR(8), Y     ),
     *                 (CURPAR(9), CHI   ), (CURPAR(10),BETA  ),
     *                 (CURPAR(11),GCM   ), (CURPAR(12),ECM   )
      COMMON /CRRANDPA/RD,FAC,U1,U2,NSEQ,ISEED,KNOR
      DOUBLE PRECISION RD(3000),FAC,U1,U2
      INTEGER          ISEED(3,10),NSEQ
      LOGICAL          KNOR

      COMMON /CRREST/  CONTNE,TAR,LT
      DOUBLE PRECISION CONTNE(3),TAR
      INTEGER          LT

      COMMON /CRRUNPAR/FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,
     *                 STEPFC,
# 4137 "corsika.h"
     *                 NRRUN,NSHOW,MPATAP,MONIIN,
     *                 MONIOU,MDEBUG,NUCNUC,MTABOUT,MLONGOUT,
     *                 ISEED1I,
# 4158 "corsika.h"
     *                 LSTCK,

     *                 LSTCK1,LSTCK2,

     *                 ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     *                 N1STTR,MDBASE,

     *                 DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     *                 FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     *                 FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN

     *                 ,FOUTFILE,IFINAM
# 4199 "corsika.h"
      COMMON /CRRUNPAC/DATDIR,DSN,DSNTAB,DSNLONG,HOST,USER
# 4215 "corsika.h"
     *                 ,FILOUT
# 4226 "corsika.h"
      DOUBLE PRECISION FIXHEI,THICK0,HILOECM,HILOELB,SIG1I,TARG1I,STEPFC

      INTEGER          NRRUN,NSHOW,MPATAP,MONIIN,MONIOU,MDEBUG,NUCNUC,
     *                 ISHOWNO,ISHW,NOPART,NRECS,NBLKS,MAXPRT,NDEBDL,
     *                 N1STTR,MDBASE,MTABOUT,MLONGOUT,ISEED1I(3)
# 4257 "corsika.h"
      INTEGER          LSTCK

     *                ,LSTCK1,LSTCK2
# 4268 "corsika.h"
      CHARACTER*132    FILOUT

      CHARACTER*255    DSN,DSNTAB,DSNLONG
      CHARACTER*132    DATDIR
      CHARACTER*20     HOST,USER
# 4288 "corsika.h"
      LOGICAL          DEBDEL,DEBUG,FDECAY,FEGS,FIRSTI,FIXINC,FIXTAR,
     *                 FIX1I,FMUADD,FNKG,FPRINT,FDBASE,FPAROUT,FTABOUT,
     *                 FLONGOUT,GHEISH,GHESIG,GHEISDB,USELOW,TMARGIN
# 4301 "corsika.h"
      LOGICAL          FOUTFILE
      INTEGER          IFINAM
# 4323 "corsika.h"

# 4333 "corsika.h"

      COMMON /CRSIGM/  SIGMA,SIGANN,SIGAIR,FRACTN,FRCTNO,SIGAIRS
      DOUBLE PRECISION SIGMA,SIGANN,SIGAIR,FRACTN,FRCTNO,SIGAIRS


      integer mmry,mxptl
      parameter (mmry=1)   !memory saving factor

C--------------------------- EPOS COMMON ------------------------------
      integer nevt,kolevt,npjevt,minfra,maxfra,nglevt,koievt,kohevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
      real phievt,bimevt,pmxevt,egyevt,xbjevt,qsqevt,zppevt,zptevt
      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra,kohevt
      parameter (mxptl=200000/mmry) !max nr of particles in nexus particle list
      integer      infragm
      common/nucl6/infragm
# 61256 "corsika.F"
      integer nptl,iorptl,idptl,ifrptl,jorptl,istptl,ibptl,ityptl
      real pptl,tivptl,xorptl
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl),jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),ityptl(mxptl)
      integer laproj,maproj,latarg,matarg
      real core,fctrmx
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
C-----------------------------------------------------------------------
      DOUBLE PRECISION EA,ELABN,ELASTI,EMAX,COSTET,PL2,PT2,PTM
CC    DOUBLE PRECISION GAMMAX
      DOUBLE PRECISION PFRX(60),PFRY(60),PTOT,CPHIV,SPHIV
      DOUBLE PRECISION FAC1,FAC2
      REAL ETOT,GNU
      INTEGER          ITYP(60),NREST,IREST,NNEW,INEW,NZNEW,NNNEW,I,
     *                 KODNEX,KODCRS,JFIN,KNEW,J,MEL,MEN,LL

      SAVE
# 61289 "corsika.F"
C-----------------------------------------------------------------------

      IF ( DEBUG ) WRITE(MDEBUG,*) 'NSTORE:'

C  NUMBER OF SPECTATORS OF REMAINING NUCLEUS IS NREST
      NREST = ITYPE/100 - npjevt
      IREST = ITYPE

      NNEW   = 0
      INEW   = 0
      ETOT   = 0.
      ELASTI = 0.
      NZNEW  = 0
      NNNEW  = 0
      KNEW   = 0

c  event variables:
c     nrevt.......... event number
c     nptevt ........ number of (stored!) particles per event
c     bimevt ........ absolute value of impact parameter
c     phievt ........ angle of impact parameter
c     kolevt ........ number of collisions
c     pmxevt ........ reference momentum
c     egyevt ........ pp cm energy (hadron) or string energy (lepton)
c     npjevt ........ number of primary projectile participants
c     ntgevt ........ number of primary target participants

      GNU  = kolevt
      EMAX = 0.D0
CC    GAMMAX = 0.D0
# 61335 "corsika.F"

C  PARTICLE LOOP
      DO  5  I = 1, nptl
C  SKIP PARTICLES NOT FROM LAST GENERATION (ISTPTL=0)
        IF ( istptl(I) .GT. 0 ) GOTO 5

c  particle variables:
c     i ............. particle number
c     idptl(i) ...... particle id
c     pptl(1,i) ..... x-component of particle momentum
c     pptl(2,i) ..... y-component of particle momentum
c     pptl(3,i) ..... z-component of particle momentum
c     pptl(4,i) ..... particle energy
c     pptl(5,i) ..... particle mass
c     iorptl(i) ..... particle number of father
c     jorptl(i) ..... particle number of mother
c     ifrptl(1,i) ... particle number of first child
c     ifrptl(2,i) ... particle number of last child
c     istptl(i) ..... generation flag: last gen.(0) or not(1) or ghost(2)
c     xorptl(1,i) ... x-component of formation point
c     xorptl(2,i) ... y-component of formation point
c     xorptl(3,i) ... z-component of formation point
c     xorptl(4,i) ... formation time
c     tivptl(1,i) ... formation time (always in the pp-cms!)
c     tivptl(2,i) ... destruction time (always in the pp-cms!)
c     ityptl(i)  .....from target (10-19), soft (20-29), hard (30-39),
c                     projectile (40-49) string, droplet (50)
c     idiptl(i) ..... id of father (999 if no father)
c     idjptl(i) ..... id of mother (999 if no mother)

C  ELIMINATE TARGET SPECTATORS
        IF ( PPTL(3,I) .LT. 0.1D0 ) GOTO 5

C  ELIMINATE BACKWARD GOING PARTICLES
        IF ( .NOT. LLONGI  .AND.  PPTL(3,I) .LT. 0. ) GOTO 5

C  CONVERT PARTICLE CODE  NEX(US) ---> C(O)RS(IKA)
C  MOST FREQUENT PARTICLES COME FIRST
        KODNEX = idptl(I)
C  NUCLEUS
        IF     ( KODNEX .GT. 10000 ) THEN
          KODCRS = MOD(KODNEX,10000)*10 + MOD(KODNEX/10000,1000)
        ELSEIF ( KODNEX .EQ.   17  ) THEN
          KODCRS = 201
        ELSEIF ( KODNEX .EQ.   18  ) THEN
          KODCRS = 301
        ELSEIF ( KODNEX .EQ.   19  ) THEN
          KODCRS = 402
C  MESONS
        ELSEIF     ( KODNEX .EQ.   110 ) THEN
          KODCRS = 7
        ELSEIF ( KODNEX .EQ.   120 ) THEN
          KODCRS = 8
        ELSEIF ( KODNEX .EQ.  -120 ) THEN
          KODCRS = 9
        ELSEIF ( KODNEX .EQ.   220 ) THEN
          KODCRS = 17
C  NUCLEONS
        ELSEIF ( KODNEX .EQ.  1220 ) THEN
          KODCRS = 13
        ELSEIF ( KODNEX .EQ.  1120 ) THEN
          KODCRS = 14
        ELSEIF ( KODNEX .EQ. -1120 ) THEN
          KODCRS = 15
        ELSEIF ( KODNEX .EQ. -1220 ) THEN
          KODCRS = 25
C  STRANGE MESONS
        ELSEIF ( KODNEX .EQ.   -20 .or. KODNEX .EQ.   230 ) THEN
          KODCRS = 10
        ELSEIF ( KODNEX .EQ.   130 ) THEN
          KODCRS = 11
        ELSEIF ( KODNEX .EQ.  -130 ) THEN
          KODCRS = 12
        ELSEIF ( KODNEX .EQ.    20 .or. KODNEX .EQ.   -230 ) THEN
          KODCRS = 16
C  STRANGE BARYONS
        ELSEIF ( KODNEX .EQ.  2130 ) THEN
          KODCRS = 18
        ELSEIF ( KODNEX .EQ.  1130 ) THEN
          KODCRS = 19
        ELSEIF ( KODNEX .EQ.  1230 ) THEN
          KODCRS = 20
        ELSEIF ( KODNEX .EQ.  2230 ) THEN
          KODCRS = 21
        ELSEIF ( KODNEX .EQ.  1330 ) THEN
          KODCRS = 22
        ELSEIF ( KODNEX .EQ.  2330 ) THEN
          KODCRS = 23
        ELSEIF ( KODNEX .EQ.  3331 ) THEN
          KODCRS = 24
        ELSEIF ( KODNEX .EQ. -2130 ) THEN
          KODCRS = 26
        ELSEIF ( KODNEX .EQ. -1130 ) THEN
          KODCRS = 27
        ELSEIF ( KODNEX .EQ. -1230 ) THEN
          KODCRS = 28
        ELSEIF ( KODNEX .EQ. -2230 ) THEN
          KODCRS = 29
        ELSEIF ( KODNEX .EQ. -1330 ) THEN
          KODCRS = 30
        ELSEIF ( KODNEX .EQ. -2330 ) THEN
          KODCRS = 31
        ELSEIF ( KODNEX .EQ. -3331 ) THEN
          KODCRS = 32
C  LEPTONS
        ELSEIF ( KODNEX .EQ.    10 ) THEN
          KODCRS = 1
        ELSEIF ( KODNEX .EQ.   -12 ) THEN
          KODCRS = 2
        ELSEIF ( KODNEX .EQ.    12 ) THEN
          KODCRS = 3
        ELSEIF ( KODNEX .EQ.   -14 ) THEN
          KODCRS = 5
        ELSEIF ( KODNEX .EQ.    14 ) THEN
          KODCRS = 6

C  CANNOT TREAT TAU LEPTONS, CONVERT TO MUONS
        ELSEIF ( KODNEX .EQ.   -16 ) THEN  ! TAU(+)
          KODCRS = 5
        ELSEIF ( KODNEX .EQ.    16 ) THEN  ! TAU(-
          KODCRS = 6
# 61544 "corsika.F"
C  CHARMED MESONS CANNOT BE TREATED, TAKE INSTEAD STRANGE MESONS
        ELSEIF ( KODNEX .EQ.  -140 ) THEN  ! D(0)
          KODCRS = 12
        ELSEIF ( KODNEX .EQ.  -240  .OR.   ! D(+)
     *           KODNEX .EQ.   240 ) THEN  ! A-D(-)
          CALL RMMARD( RD,1,1 )
          IF ( RD(1) .GE. 0.5D0 ) THEN
            KODCRS = 10
          ELSE
            KODCRS = 16
          ENDIF
        ELSEIF ( KODNEX .EQ.   140 ) THEN  ! A-D(0)
          KODCRS = 11
        ELSEIF ( KODNEX .EQ.  -340  .OR.   ! D_S(+)
     *           KODNEX .EQ.   340 ) THEN  ! A-D_S(-)
          IF ( RD(1) .GE. 0.5D0 ) THEN
            KODCRS = 10
          ELSE
            KODCRS = 16
          ENDIF
        ELSEIF ( KODNEX .EQ.   440 ) THEN  ! ETA_C
          KODCRS = 17
        ELSEIF ( KODNEX .EQ.  -141 ) THEN  ! D*(0)
          KODCRS = 12
        ELSEIF ( KODNEX .EQ.  -241  .OR.   ! D*(+)
     *           KODNEX .EQ.   241  .OR.   ! A-D*(-)
     *           KODNEX .EQ.  -341  .OR.   ! D_S*(+)
     *           KODNEX .EQ.   341 ) THEN  ! A-D_S*(-)
          IF ( RD(1) .GE. 0.5D0 ) THEN
            KODCRS = 10
          ELSE
            KODCRS = 16
          ENDIF
        ELSEIF ( KODNEX .EQ.   141 ) THEN  ! A-D*(0)
          KODCRS = 11
        ELSEIF ( KODNEX .EQ.   441 ) THEN  ! J/PSI
          KODCRS = 17
C  CHARMED BARYONS CANNOT BE TREATED, TAKE INSTEAD STRANGE BARYON
        ELSEIF ( KODNEX .EQ.  2140 ) THEN  ! LAMBDA_C(+)
          KODCRS = 18
        ELSEIF ( KODNEX .EQ.  3140 ) THEN  ! XI_C(+)
          KODCRS = 22
        ELSEIF ( KODNEX .EQ.  3240 ) THEN  ! XI_C(0)
          KODCRS = 23
        ELSEIF ( KODNEX .EQ.  1140 ) THEN  ! SIGMA_C(++)
          KODCRS = 19
        ELSEIF ( KODNEX .EQ.  1240 ) THEN  ! SIGMA_C(+)
          KODCRS = 20
        ELSEIF ( KODNEX .EQ.  2240 ) THEN  ! SIGMA_C(0)
          KODCRS = 21
        ELSEIF ( KODNEX .EQ.  1340 ) THEN  ! XI_C''(+)
          KODCRS = 23
        ELSEIF ( KODNEX .EQ.  2340 ) THEN  ! XI_C''(0)
          KODCRS = 22
        ELSEIF ( KODNEX .EQ.  3340 ) THEN  ! OMEGA_C(0)
          KODCRS = 24
        ELSEIF ( KODNEX .EQ. -2140 ) THEN  ! A-LAMBDA_C(-)
          KODCRS = 26
        ELSEIF ( KODNEX .EQ. -3140 ) THEN  ! A-XI_C(-)
          KODCRS = 30
        ELSEIF ( KODNEX .EQ. -3240 ) THEN  ! A-XI_C(0)
          KODCRS = 31
        ELSEIF ( KODNEX .EQ. -1140 ) THEN  ! A-SIGMA_C(--)
          KODCRS = 27
        ELSEIF ( KODNEX .EQ. -1240 ) THEN  ! A-SIGMA_C(-)
          KODCRS = 28
        ELSEIF ( KODNEX .EQ. -2240 ) THEN  ! A-SIGMA_C(0)
          KODCRS = 29
        ELSEIF ( KODNEX .EQ. -1340 ) THEN  ! A-XI_C''(-)
          KODCRS = 31
        ELSEIF ( KODNEX .EQ. -2340 ) THEN  ! A-XI_C''(0)
          KODCRS = 30
        ELSEIF ( KODNEX .EQ. -3340 ) THEN  ! A-OMEGA_C(0)
          KODCRS = 32
        ELSEIF ( KODNEX .EQ.  1141 ) THEN  ! SIGMA_C*(++)
          KODCRS = 161
        ELSEIF ( KODNEX .EQ.  1241 ) THEN  ! SIGMA_C*(+)
          KODCRS = 162
        ELSEIF ( KODNEX .EQ.  2241 ) THEN  ! SIGMA_C*(0)
          KODCRS = 163
        ELSEIF ( KODNEX .EQ. -1141 ) THEN  ! A-SIGMA_C*(--)
          KODCRS = 171
        ELSEIF ( KODNEX .EQ. -1241 ) THEN  ! A-SIGMA_C*(-)
          KODCRS = 172
        ELSEIF ( KODNEX .EQ. -2241 ) THEN  ! A-SIGMA_C*(0)
          KODCRS = 173
# 61654 "corsika.F"
C  NEUTRINOS ARE SKIPPED
        ELSEIF ( KODNEX .EQ.    11 ) THEN
          GOTO 55
        ELSEIF ( KODNEX .EQ.   -11 ) THEN
          GOTO 55
        ELSEIF ( KODNEX .EQ.    13 ) THEN
          GOTO 55
        ELSEIF ( KODNEX .EQ.   -13 ) THEN
          GOTO 55
        ELSEIF ( KODNEX .EQ.    15 ) THEN
          GOTO 55
        ELSEIF ( KODNEX .EQ.   -15 ) THEN
          GOTO 55

        ELSE
          WRITE(MONIOU,*)'NSTORE: UNKNOWN PARTICLE CODE IDPTL=',idptl(I)
          GOTO 5
        ENDIF
        SECPAR(0) = KODCRS



      RETURN
      END

*-- Author :    D. HECK IK FZK KARLSRUHE   18/03/2003
C=======================================================================

      FUNCTION RANGEN()

C-----------------------------------------------------------------------
C  RAN(DOM  NUMBER) GEN(ERATOR)
C
C  SEE SUBROUT. RMMARD
C  WE USE HERE A SIMPLIFIED FORM OF RMMARD WITH JSEQ=1, LENV=1.
C  THIS FUNCTION IS CALLED FROM MANY EPOS/NEXUS ROUTINES.
C-----------------------------------------------------------------------

      IMPLICIT NONE


      INTEGER          KSEQ
      PARAMETER        (KSEQ = 8)
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS

      COMMON /CRRANMA4/C,U,IJKL,I97,J97,NTOT,NTOT2,JSEQ
      DOUBLE PRECISION C(KSEQ),U(97,KSEQ),UNI
      INTEGER          IJKL(KSEQ),I97(KSEQ),J97(KSEQ),
     *                 NTOT(KSEQ),NTOT2(KSEQ),JSEQ

    

      COMMON /CXINPUT/ CXTHR,CXMCC,CXWMX,CXMCS

     *                ,INLUN,OUTLUN,IDCX,ISX0,FINCNX,FCXGHE,FCXWMX,FCXCE
     *                ,FCORS
      DOUBLE PRECISION CXTHR(3),CXMCC(3),CXWMX(3),CXMCS

      INTEGER          INLUN,OUTLUN,IDCX,ISX0
      LOGICAL          FINCNX,FCXGHE,FCXWMX,FCXCE,FCORS
      COMMON /CXCONVE/ CXXCONV,CXYCONV,CXTCONV
      DOUBLE PRECISION CXXCONV,CXYCONV,CXTCONV

      REAL             RANGEN
      SAVE
C-----------------------------------------------------------------------

      JSEQ = 1

 1    CONTINUE
      UNI = U(I97(JSEQ),JSEQ) - U(J97(JSEQ),JSEQ)
      IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
      U(I97(JSEQ),JSEQ) = UNI
      I97(JSEQ)  = I97(JSEQ) - 1
      IF ( I97(JSEQ) .EQ. 0 ) I97(JSEQ) = 97
      J97(JSEQ)  = J97(JSEQ) - 1
      IF ( J97(JSEQ) .EQ. 0 ) J97(JSEQ) = 97
      C(JSEQ)    = C(JSEQ) - CD
      IF ( C(JSEQ) .LT. 0.D0 ) C(JSEQ)  = C(JSEQ) + CM
      UNI        = UNI - C(JSEQ)
      IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
C  AN EXACT ZERO HERE IS VERY UNLIKELY, BUT LET''S BE SAFE.
      IF ( UNI .EQ. 0.D0 ) UNI = TWOM48
      RANGEN = SNGL( UNI )

      NTOT(JSEQ) = NTOT(JSEQ) + 1
      IF ( NTOT(JSEQ) .GE. MODCNS )  THEN
        NTOT2(JSEQ) = NTOT2(JSEQ) + 1
        NTOT(JSEQ)  = NTOT(JSEQ) - MODCNS
      ENDIF
C  AN EXACT ZERO HERE IS VERY UNLIKELY, BUT LET''S BE SAFE AND
C  TAKE A NEW RANDOM NUMBER
      IF     ( RANGEN .EQ. 0. ) THEN
        GO TO 1
      ELSEIF ( RANGEN .EQ. 1. ) THEN
        GO TO 1
      ENDIF

      RETURN
      END

*-- Author :    T. PIEROG IK FZK KARLSRUHE   24/11/2005
C=======================================================================

      DOUBLE PRECISION FUNCTION G900GT()

C-----------------------------------------------------------------------
C  GET SEED FROM RANDOM  NUMBER GENERATOR
C
C  THIS FUNCTION IS CALLED FROM RANFGT IN EPOS/NEXUS ROUTINES.
C-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER          KSEQ
      PARAMETER        (KSEQ = 8)
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS

      COMMON /CRRANMA4/C,U,IJKL,I97,J97,NTOT,NTOT2,JSEQ
      DOUBLE PRECISION C(KSEQ),U(97,KSEQ),UNI
      INTEGER          IJKL(KSEQ),I97(KSEQ),J97(KSEQ),
     *                 NTOT(KSEQ),NTOT2(KSEQ),JSEQ


      COMMON /CXINPUT/ CXTHR,CXMCC,CXWMX,CXMCS

     *                ,INLUN,OUTLUN,IDCX,ISX0,FINCNX,FCXGHE,FCXWMX,FCXCE
     *                ,FCORS
      DOUBLE PRECISION CXTHR(3),CXMCC(3),CXWMX(3),CXMCS

      INTEGER          INLUN,OUTLUN,IDCX,ISX0
      LOGICAL          FINCNX,FCXGHE,FCXWMX,FCXCE,FCORS
      COMMON /CXCONVE/ CXXCONV,CXYCONV,CXTCONV
      DOUBLE PRECISION CXXCONV,CXYCONV,CXTCONV

     

      SAVE
C-----------------------------------------------------------------------
      JSEQ = 1
# 62332 "corsika.F"
      UNI = NTOT(JSEQ) + 1.D9*NTOT2(JSEQ)
      G900GT = UNI

      RETURN
      END

*-- Author :    T. PIEROG IK FZK KARLSRUHE   24/11/2005
C=======================================================================

      DOUBLE PRECISION FUNCTION G900ST(DSEED)

C-----------------------------------------------------------------------
C  STORE SEED FROM RANDOM  NUMBER GENERATOR (DUMMY FUNCTION IN CORSIKA)
C
C  THIS FUNCTION IS CALLED FROM RANFST IN EPOS/NEXUS ROUTINES.
C-----------------------------------------------------------------------

      IMPLICIT NONE


      INTEGER          KSEQ
      PARAMETER        (KSEQ = 8)
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS

      COMMON /CRRANMA4/C,U,IJKL,I97,J97,NTOT,NTOT2,JSEQ
      DOUBLE PRECISION C(KSEQ),U(97,KSEQ),UNI
      INTEGER          IJKL(KSEQ),I97(KSEQ),J97(KSEQ),
     *                 NTOT(KSEQ),NTOT2(KSEQ),JSEQ

     

      COMMON /CXINPUT/ CXTHR,CXMCC,CXWMX,CXMCS

     *                ,INLUN,OUTLUN,IDCX,ISX0,FINCNX,FCXGHE,FCXWMX,FCXCE
     *                ,FCORS
      DOUBLE PRECISION CXTHR(3),CXMCC(3),CXWMX(3),CXMCS

      INTEGER          INLUN,OUTLUN,IDCX,ISX0
      LOGICAL          FINCNX,FCXGHE,FCXWMX,FCXCE,FCORS
      COMMON /CXCONVE/ CXXCONV,CXYCONV,CXTCONV
      DOUBLE PRECISION CXXCONV,CXYCONV,CXTCONV

     
# 62354 "corsika.F" 2

      DOUBLE PRECISION DSEED
      SAVE
C-----------------------------------------------------------------------
      JSEQ = 1
# 62378 "corsika.F"
      UNI    = DSEED
      G900ST = UNI

      RETURN
      END

*-- Author :    T. PIEROG IK FZK KARLSRUHE   24/11/2005
C=======================================================================

      SUBROUTINE RANFGT(SEED)

C-----------------------------------------------------------------------

      DOUBLE PRECISION SEED,G900GT,G900ST,DUMMY
C-----------------------------------------------------------------------
      SEED  = G900GT()
      RETURN

      ENTRY RANFST(SEED)
      DUMMY = G900ST(SEED)

      RETURN
      END

*-- Author :    D. HECK IK FZK KARLSRUHE   18/03/2003
C=======================================================================

      DOUBLE PRECISION FUNCTION DRANGEN(DUMMY)

C-----------------------------------------------------------------------
C  D(OUBLE PRECISION) RAN(DOM  NUMBER) GEN(ERATOR)
C
C  SEE SUBROUT. RMMARD
C  WE USE HERE A SIMPLIFIED FORM OF RMMARD WITH JSEQ=1, LENV=1.
C  THIS FUNCTION IS CALLED FROM MANY EPOS ROUTINES.
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER          KSEQ
      PARAMETER        (KSEQ = 8)
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS

      COMMON /CRRANMA4/C,U,IJKL,I97,J97,NTOT,NTOT2,JSEQ
      DOUBLE PRECISION C(KSEQ),U(97,KSEQ),UNI
      INTEGER          IJKL(KSEQ),I97(KSEQ),J97(KSEQ),
     *                 NTOT(KSEQ),NTOT2(KSEQ),JSEQ


      COMMON /CXINPUT/ CXTHR,CXMCC,CXWMX,CXMCS

     *                ,INLUN,OUTLUN,IDCX,ISX0,FINCNX,FCXGHE,FCXWMX,FCXCE
     *                ,FCORS
      DOUBLE PRECISION CXTHR(3),CXMCC(3),CXWMX(3),CXMCS

      INTEGER          INLUN,OUTLUN,IDCX,ISX0
      LOGICAL          FINCNX,FCXGHE,FCXWMX,FCXCE,FCORS
      COMMON /CXCONVE/ CXXCONV,CXYCONV,CXTCONV
      DOUBLE PRECISION CXXCONV,CXYCONV,CXTCONV


      DOUBLE PRECISION DUMMY
      SAVE
C-----------------------------------------------------------------------

      JSEQ = 1

 1    CONTINUE
      UNI = U(I97(JSEQ),JSEQ) - U(J97(JSEQ),JSEQ)
      IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
      U(I97(JSEQ),JSEQ) = UNI
      I97(JSEQ)  = I97(JSEQ) - 1
      IF ( I97(JSEQ) .EQ. 0 ) I97(JSEQ) = 97
      J97(JSEQ)  = J97(JSEQ) - 1
      IF ( J97(JSEQ) .EQ. 0 ) J97(JSEQ) = 97
      C(JSEQ)    = C(JSEQ) - CD
      IF ( C(JSEQ) .LT. 0.D0 ) C(JSEQ)  = C(JSEQ) + CM
      UNI        = UNI - C(JSEQ)
      IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
C  AN EXACT ZERO HERE IS VERY UNLIKELY, BUT LET''S BE SAFE.
      IF ( UNI .EQ. 0.D0 ) UNI = TWOM48
      DRANGEN = UNI

      NTOT(JSEQ) = NTOT(JSEQ) + 1
      IF ( NTOT(JSEQ) .GE. MODCNS )  THEN
        NTOT2(JSEQ) = NTOT2(JSEQ) + 1
        NTOT(JSEQ)  = NTOT(JSEQ) - MODCNS
      ENDIF

      RETURN
      END

# 62514 "corsika.F"

C=======================================================================

      SUBROUTINE RANFCV(SEED)

C-----------------------------------------------------------------------
c  Dummy function in CORSIKA
c  Convert input seed to EPOS random number seed
c  Since input seed and EPOS (from Corsika) seed are different,
c  define input seed as : seed=ISEED(3)*1E9+ISEED(2)
c-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION SEED,DUMMY

      DUMMY = SEED

      RETURN
      END

C=======================================================================

      SUBROUTINE RANFINI(SEED,ISEQ,IQQ)

c-----------------------------------------------------------------------
c  Dummy function in CORSIKA
c  Initialize random number sequence iseq with seed
c  if iqq=-1, run first ini
c    iqq=0 , set what sequence should be used
c    iqq=1 , initialize sequence for initialization
c    iqq=2 , initialize sequence for first event
c-----------------------------------------------------------------------

      IMPLICIT NONE
      DOUBLE PRECISION SEED,DUMMY
      INTEGER          ISEQ,IQQ,IDUM

      DUMMY = SEED
      IDUM  = IQQ
      IDUM  = ISEQ

      RETURN
      END

c-----------------------------------------------------------------------
      subroutine urqmd(idum)
c-----------------------------------------------------------------------
c Dummy function for compilation of EPOS if URQMD is not selected
c Initialize random number sequence iseq with seed
c-----------------------------------------------------------------------

      IMPLICIT NONE
      integer idum0,idum
 
      idum0=idum

      return
      end

