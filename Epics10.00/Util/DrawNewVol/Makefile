include  $(EPICSTOP)/site.config


EXTHDRS	      = ../../prog/NewVol/Zciecone.h \
		$(COSMOSINC)/Zglobalc.h \
		$(EPICSINC)/Zcnfig.h \
		$(EPICSINC)/Zcnfig2.h \
		$(EPICSINC)/Zep3Vec.h \
		$(EPICSINC)/ZepDirec.h \
		$(EPICSINC)/ZepDraw.h \
		$(EPICSINC)/ZepMaxdef.h \
		$(EPICSINC)/ZepPos.h

HDRS	      = Zelightg.h


MAKEFILE      = Makefile

OBJS	      = epDraw_angle.o \
		epDraw_ccone.o \
		epDraw_ccyl.o \
		epDraw_cecyl.o \
		epDraw_chute.o \
		epDraw_ciecone.o \
		epDraw_cpipe.o \
		epDraw_csph.o \
		epDraw_dum1.o \
		epDraw_dum2.o \
		epDraw_dum3.o \
		epDraw_dum4.o \
		epDraw_dum5.o \
		epDraw_dum6.o \
		epDraw_dum7.o \
		epDraw_dum8.o \
		epDraw_dum9.o \
		epDraw_elightg.o \
		epDraw_ellips.o \
		epDraw_eparab.o \
		epDraw_gbox.o \
		epDraw_horse.o \
		epDraw_hypb.o \
		epDraw_hypb2.o \
		epDraw_ocone.o \
		epDraw_ocyl.o \
		epDraw_oecyl.o \
		epDraw_opipe.o \
		epDraw_sccyl.o \
		epDraw_scyl.o \
		epDraw_sqpipe.o \
		epDraw_swave.o \
		epDraw_tripyra.o \
		epDraw_wave.o \
		epDraw_octagon.o \
		epDraw_honeycomb.o  \
		epDraw_sqTccl.o  \
		epDraw_torus.o  \
		epDraw_fpolygon.o \
		epDraw_polyhed0.o
PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = epDraw_angle.f \
		epDraw_ccone.f \
		epDraw_ccyl.f \
		epDraw_cecyl.f \
		epDraw_chute.f \
		epDraw_ciecone.f \
		epDraw_cpipe.f \
		epDraw_csph.f \
		epDraw_dum1.f \
		epDraw_dum2.f \
		epDraw_dum3.f \
		epDraw_dum4.f \
		epDraw_dum5.f \
		epDraw_dum6.f \
		epDraw_dum7.f \
		epDraw_dum8.f \
		epDraw_dum9.f \
		epDraw_elightg.f \
		epDraw_ellips.f \
		epDraw_eparab.f \
		epDraw_gbox.f \
		epDraw_horse.f \
		epDraw_hypb.f \
		epDraw_hypb2.f \
		epDraw_ocone.f \
		epDraw_ocyl.f \
		epDraw_oecyl.f \
		epDraw_opipe.f \
		epDraw_sccyl.f \
		epDraw_scyl.f \
		epDraw_sqpipe.f \
		epDraw_swave.f \
		epDraw_tripyra.f \
		epDraw_wave.f  \
		epDraw_octagon.f \
		epDraw_honeycomb.f \
		epDraw_sqTccl.f \
		epDraw_torus.f  \
		epDraw_fpolygon.f \
		epDraw_polyhed0.f 
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
epDraw_angle.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_ccone.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_ccyl.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_cecyl.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_chute.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_ciecone.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h ../../prog/NewVol/Zciecone.h
epDraw_cpipe.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h \
	$(COSMOSINC)/Zglobalc.h
epDraw_csph.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_dum1.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_dum2.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_dum3.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_dum4.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_dum5.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_dum6.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_dum7.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_dum8.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_dum9.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_elightg.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h Zelightg.h
epDraw_ellips.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_eparab.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_gbox.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_horse.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_hypb.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_hypb2.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_ocone.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_ocyl.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_oecyl.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_opipe.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_sccyl.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h \
	$(EPICSINC)/ZepDirec.h
epDraw_scyl.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_sqpipe.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_swave.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_tripyra.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_wave.o: $(COSMOSINC)/Zglobalc.h \
	$(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_octagon.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h
epDraw_honeycomb.o: $(EPICSINC)/Zep3Vec.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/ZepMaxdef.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/ZepDraw.h \
	$(EPICSINC)/ZepPos.h

$(LIBRARY):         \
               $(LIBRARY)(epDraw_angle.o) \
               $(LIBRARY)(epDraw_ccone.o) \
               $(LIBRARY)(epDraw_ccyl.o) \
               $(LIBRARY)(epDraw_cecyl.o) \
               $(LIBRARY)(epDraw_chute.o) \
               $(LIBRARY)(epDraw_ciecone.o) \
               $(LIBRARY)(epDraw_cpipe.o) \
               $(LIBRARY)(epDraw_csph.o) \
               $(LIBRARY)(epDraw_dum1.o) \
               $(LIBRARY)(epDraw_dum2.o) \
               $(LIBRARY)(epDraw_dum3.o) \
               $(LIBRARY)(epDraw_dum4.o) \
               $(LIBRARY)(epDraw_dum5.o) \
               $(LIBRARY)(epDraw_dum6.o) \
               $(LIBRARY)(epDraw_dum7.o) \
               $(LIBRARY)(epDraw_dum8.o) \
               $(LIBRARY)(epDraw_dum9.o) \
               $(LIBRARY)(epDraw_elightg.o) \
               $(LIBRARY)(epDraw_ellips.o) \
               $(LIBRARY)(epDraw_eparab.o) \
               $(LIBRARY)(epDraw_gbox.o) \
               $(LIBRARY)(epDraw_horse.o) \
               $(LIBRARY)(epDraw_hypb.o) \
               $(LIBRARY)(epDraw_hypb2.o) \
               $(LIBRARY)(epDraw_ocone.o) \
               $(LIBRARY)(epDraw_ocyl.o) \
               $(LIBRARY)(epDraw_oecyl.o) \
               $(LIBRARY)(epDraw_opipe.o) \
               $(LIBRARY)(epDraw_sccyl.o) \
               $(LIBRARY)(epDraw_scyl.o) \
               $(LIBRARY)(epDraw_sqpipe.o) \
               $(LIBRARY)(epDraw_swave.o) \
               $(LIBRARY)(epDraw_tripyra.o) \
               $(LIBRARY)(epDraw_wave.o)  \
               $(LIBRARY)(epDraw_octagon.o) \
               $(LIBRARY)(epDraw_honeycomb.o) \
               $(LIBRARY)(epDraw_sqTccl.o) \
               $(LIBRARY)(epDraw_torus.o) \
               $(LIBRARY)(epDraw_fpolygon.o) \
               $(LIBRARY)(epDraw_polyhed0.o)
