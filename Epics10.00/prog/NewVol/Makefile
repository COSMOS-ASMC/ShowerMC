include  $(EPICSTOP)/site.config


EXTHDRS	      = $(COSMOSINC)/Zglobalc.h \
		$(EPICSINC)/Zcnfig.h \
		$(EPICSINC)/Zcnfig2.h \
		$(EPICSINC)/Zep3Vec.h \
		$(EPICSINC)/ZepDirec.h \
		$(EPICSINC)/ZepMaxdef.h \
		$(EPICSINC)/ZepPos.h \
		$(EPICSINC)/ZepTrackp.h \
		$(EPICSINC)/Zepdebug.h

HDRS	      = Zciecone.h \
		Zelightg.h


MAKEFILE      = Makefile

OBJS	      = ep_angle.o \
		ep_ccone.o \
		ep_ccyl.o \
		ep_cecyl.o \
		ep_chute.o \
		ep_ciecone.o \
		ep_cpipe.o \
		ep_csph.o \
		ep_dcsph.o \
		ep_elightg.o \
		ep_ellips.o \
		ep_eparab.o \
		ep_gbox.o \
		ep_horse.o \
		ep_hypb.o \
		ep_hypb2.o \
		ep_ocone.o \
		ep_ocyl.o \
		ep_oecyl.o \
		ep_opipe.o \
		ep_sccyl.o \
		ep_scyl.o \
		ep_sqpipe.o \
		ep_swave.o \
		ep_tripyra.o \
		ep_wave.o \
		ep_honeycomb.o \
		ep_sqtccl.o \
		epdummy1.o \
		epdummy2.o \
		epdummy3.o \
		epdummy4.o \
		epdummy5.o \
		epdummy6.o \
		epdummy7.o \
		epdummy8.o \
		epdummy9.o \
		epxpcylwall.o \
		epxpellipse.o  \
		ep_torus.o \
		ep_fpolygon.o \
		ep_polyhed0.o
PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = ep_angle.f \
		ep_ccone.f \
		ep_ccyl.f \
		ep_cecyl.f \
		ep_chute.f \
		ep_ciecone.f \
		ep_cpipe.f \
		ep_csph.f \
		ep_dcsph.f \
		ep_elightg.f \
		ep_ellips.f \
		ep_eparab.f \
		ep_gbox.f \
		ep_horse.f \
		ep_hypb.f \
		ep_hypb2.f \
		ep_ocone.f \
		ep_ocyl.f \
		ep_oecyl.f \
		ep_opipe.f \
		ep_sccyl.f \
		ep_scyl.f \
		ep_sqpipe.f \
		ep_swave.f \
		ep_tripyra.f \
		ep_wave.f \
		ep_honeycomb.f \
		ep_sqtccl.f \
		epdummy1.f \
		epdummy2.f \
		epdummy3.f \
		epdummy4.f \
		epdummy5.f \
		epdummy6.f \
		epdummy7.f \
		epdummy8.f \
		epdummy9.f \
		epxpcylwall.f \
		epxpellipse.f \
		ep_torus.f \
		ep_fpolygon.f \
		ep_polyhed0.f 
SYSHDRS	      =

all:
	make	Checkdir $(LIBRARY) SeeSubDir "TARGET=$@"

clean:
	@rm -f $(OBJS) core  *~ a.out 
	make SeeSubDir "TARGET=$@"

veryclean:
	@rm -f $(OBJS) core  *~ a.out temp*.f
	make SeeSubDir "TARGET=$@"

clobber:
	@rm -f $(OBJS) $(LIBRARY) core tags
	make SeeSubDir "TARGET=$@"

depend:
	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)
	make SeeSubDir "TARGET=$@"

echo:
	@echo $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"

extract:
	@ar x $(LIBRARY)


index:
	@ctags -wx $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"

print:
	@$(PRINT) $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"

tags:           $(HDRS) $(SRCS)
	@ctags $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"


Checkdir:
	@( cd $(EPICSTOP)/lib; if [ ! -d  $(ARCH) ]; then  mkdir $(ARCH); fi)


SeeSubDir:
	@for i in  * \
	;do   \
		if [ -d  $$i ] ; then \
			(cd $$i; \
			if [ -f Makefile  -a  $$i != Test \
			-a $$i != Util -a $$i != UC  \
	                -a $$i != Data  -a $$i != epics \
	                -a $$i != Obso  -a $$i != TeX \
	                -a $$i != lib  -a $$i != UserHook \
	                -a $$i != Scrpt -a $$i != DistJob ] ; then \
		            $(MAKE) $(TARGET) \
			;fi) \
		; fi \
	done

###
ep_angle.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h \
	$(COSMOSINC)/Zglobalc.h
ep_ccone.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_ccyl.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/ZepPos.h
ep_cecyl.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/ZepPos.h
ep_chute.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_ciecone.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h Zciecone.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_cpipe.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
ep_csph.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/ZepPos.h
ep_dcsph.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/ZepPos.h
ep_elightg.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h Zelightg.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_ellips.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_eparab.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_gbox.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_horse.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_hypb.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_hypb2.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_ocone.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_ocyl.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_oecyl.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_opipe.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_sccyl.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h
ep_scyl.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h
ep_sqpipe.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_swave.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_tripyra.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
ep_wave.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h

ep_honeycomb.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h \
	$(EPICSINC)/Zepdebug.h
$(LIBRARY)(epdummy1.o): $(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
$(LIBRARY)(epdummy2.o): $(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
$(LIBRARY)(epdummy3.o): $(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
$(LIBRARY)(epdummy4.o): $(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
$(LIBRARY)(epdummy5.o): $(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
$(LIBRARY)(epdummy6.o): $(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
$(LIBRARY)(epdummy7.o): $(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
$(LIBRARY)(epdummy8.o): $(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
$(LIBRARY)(epdummy9.o): $(EPICSINC)/ZepTrackp.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
$(LIBRARY)(epxpcylwall.o): $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
$(LIBRARY)(ep_torus.o): $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
$(LIBRARY)(ep_fpolygon.o): $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h

$(LIBRARY)(ep_polyhed0.o): $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h

$(LIBRARY):         \
               $(LIBRARY)(ep_angle.o) \
               $(LIBRARY)(ep_ccone.o) \
               $(LIBRARY)(ep_ccyl.o) \
               $(LIBRARY)(ep_cecyl.o) \
               $(LIBRARY)(ep_chute.o) \
               $(LIBRARY)(ep_ciecone.o) \
               $(LIBRARY)(ep_cpipe.o) \
               $(LIBRARY)(ep_csph.o) \
               $(LIBRARY)(ep_dcsph.o) \
               $(LIBRARY)(ep_elightg.o) \
               $(LIBRARY)(ep_ellips.o) \
               $(LIBRARY)(ep_eparab.o) \
               $(LIBRARY)(ep_gbox.o) \
               $(LIBRARY)(ep_horse.o) \
               $(LIBRARY)(ep_hypb.o) \
               $(LIBRARY)(ep_hypb2.o) \
               $(LIBRARY)(ep_ocone.o) \
               $(LIBRARY)(ep_ocyl.o) \
               $(LIBRARY)(ep_oecyl.o) \
               $(LIBRARY)(ep_opipe.o) \
               $(LIBRARY)(ep_sccyl.o) \
               $(LIBRARY)(ep_scyl.o) \
               $(LIBRARY)(ep_sqpipe.o) \
               $(LIBRARY)(ep_swave.o) \
               $(LIBRARY)(ep_tripyra.o) \
               $(LIBRARY)(ep_wave.o) \
               $(LIBRARY)(ep_honeycomb.o) \
               $(LIBRARY)(ep_sqtccl.o) \
               $(LIBRARY)(epdummy1.o) \
               $(LIBRARY)(epdummy2.o) \
               $(LIBRARY)(epdummy3.o) \
               $(LIBRARY)(epdummy4.o) \
               $(LIBRARY)(epdummy5.o) \
               $(LIBRARY)(epdummy6.o) \
               $(LIBRARY)(epdummy7.o) \
               $(LIBRARY)(epdummy8.o) \
               $(LIBRARY)(epdummy9.o) \
               $(LIBRARY)(epxpcylwall.o) \
               $(LIBRARY)(epxpellipse.o)  \
               $(LIBRARY)(ep_torus.o) \
               $(LIBRARY)(ep_fpolygon.o) \
               $(LIBRARY)(ep_polyhed0.o)
