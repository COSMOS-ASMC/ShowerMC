!======================================================================c
!...Switches and parameters of JAM
      common/jamdat1/mstc(200),parc(200),mstd(200),pard(200)
      common/jamdat2/mste(50),pare(50)
      save /jamdat1/,/jamdat2/

!...Output file names and message.
      common/jamdat3/fname(8)
      character fname*80
      save /jamdat3/
      common/jamerdat/check(10)
      character check*120
      save /jamerdat/

!...PYTHIA Commonblocks.
      common/jydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/jydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      common/jydat3/mdcy(500,3),mdme(4000,3),brat(4000),kfdp(4000,5)
      common/jydat4/chaf(500,2)
      character chaf*16
      save /jydat1/,/jydat2/,/jydat3/,/jydat4/

!...Elementary particles.
      common/partdat3/id_quark,id_lept,id_exc,id_boson,id_diq
     $       ,id_tec,id_susy,id_special
!...Mesons.
      common/partdat2/id_pi,id_light1,id_light0
     $    ,id_str,id_charm,id_bott,id_cc,id_bb,id_mdiff
!...Baryons.
      common/partdat1/id_nucl,id_nucls,id_delt,id_delts
     $     ,id_lamb,id_lambs,id_sigm,id_sigms,id_xi,id_xis,id_omega
     $    ,id_charmb,id_bottb,id_bdiff
      save /partdat1/,/partdat2/,/partdat3/

!======================================================================c
