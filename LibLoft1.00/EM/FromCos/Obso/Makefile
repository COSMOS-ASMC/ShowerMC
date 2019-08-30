include  $(COSMOSTOP)/site.config


EXTHDRS	      = $(COSMOSINC)/ZbpCnst.h \
		$(COSMOSINC)/ZbpTable.h \
		$(COSMOSINC)/Zcode.h \
		$(COSMOSINC)/Zcondc.h \
		$(COSMOSINC)/Zcoord.h \
		$(COSMOSINC)/ZdedxAir.h \
		$(COSMOSINC)/Zdirec.h \
		$(COSMOSINC)/Zelemagp.h \
		$(COSMOSINC)/Zglobalc.h \
		$(COSMOSINC)/Zmagfield.h \
		$(COSMOSINC)/Zmass.h \
		$(COSMOSINC)/Zmaxdef.h \
		$(COSMOSINC)/Zpos.h \
		$(COSMOSINC)/Zptcl.h \
		$(COSMOSINC)/Ztrack.h \
		$(COSMOSINC)/Ztrackv.h \
		$(COSMOSINC)/Zunionmap.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = cBrCSamp.o \
		cBrLSamp.o \
		cBrSSamp.o \
		cSampMolRedA.o \
		cGetXXsec.o \
		cKnock.o \
		cPrCSamp.o \
		cPrLSamp.o \
		cPrTSamp.o \
		canihiPath.o \
		cbhabhaPath.o \
		cbremErgLPM.o \
		cbremLMPXsec.o \
		cbremsPath.o \
		ccomptPath.o \
		cdedx2atLE.o \
		cdedxElec.o \
		cdedxEleci.o \
		cdedxInAir.o \
		cdedxNone.o \
		cdedxatHE.o \
		cdedxdlt.o \
		cdedxe.o \
		cexcessLen.o \
		cmigdFunc.o \
		cmoliere.o \
		cmollerPath.o \
		cpairErgLPM.o \
		cpairLPMXsec.o \
		cpairPath.o \
		cphotoEEPath.o \
		cMolfunc.o \
		cSampMol.o 

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = cBrCSamp.f \
		cBrLSamp.f \
		cBrSSamp.f \
		cSampMolRedA.f90 \
		cGetXXsec.f \
		cKnock.f \
		cPrCSamp.f \
		cPrLSamp.f \
		cPrTSamp.f \
		canihiPath.f \
		cbhabhaPath.f \
		cbremErgLPM.f \
		cbremLMPXsec.f \
		cbremsPath.f \
		ccomptPath.f \
		cdedx2atLE.f \
		cdedxElec.f \
		cdedxEleci.f \
		cdedxInAir.f \
		cdedxNone.f \
		cdedxatHE.f \
		cdedxdlt.f \
		cdedxe.f \
		cexcessLen.f \
		cmigdFunc.f \
		cmoliere.f \
		cmollerPath.f \
		cpairErgLPM.f \
		cpairLPMXsec.f \
		cpairPath.f \
		cphotoEEPath.f \
		cMolfunc.f90 \
		cSampMol.f

SYSHDRS	      =

all:
	make	Checkdir $(LIBRARY) SeeSubDir "TARGET=$@"

clean:
	@rm -f $(OBJS) core  *~ a.out temp*.o 
	make SeeSubDir "TARGET=$@"

veryclean:
	@rm -f $(OBJS) core  *~ a.out temp*.o temp*.f
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

install:	$(LIBRARY)
		@echo "****" Installing $(LIBRARY) into $(DEST)
		$(INSTALL) -f  $(DEST) $(LIBRARY)
		@echo "****" Installing Template files into $(TEMPLATE)
		cp -r  $(COSMOSTOP)/Data  $(COSMOSTOP)/Util *
		 $(COSMOSTOP)/UserHook  $(COSMOSTOP)/Scrpt *
		 $(COSMOSTOP)/Doc        $(TEMPLATE)/
		cp $(COSMOSTOP)/site.config $(TEMPLATE)/
		@echo "****" Installing header files into $(HEADER)
		cp -r $(COSMOSTOP)/cosmos  $(HEADER)/
		@echo "****" Installation ended.
		@echo "****" Owner and group should be changed by hand

print:
	@$(PRINT) $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"

tags:           $(HDRS) $(SRCS)
	@ctags $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"


Checkdir:
	@( cd $(COSMOSTOP)/lib; if [ ! -d  $(ARCH) ]; then  mkdir $(ARCH); fi)



SeeSubDir:
	@for i in TPXS DCS * \
	;do   \
		if [ -d  $$i ] ; then \
			(cd $$i; \
			if [ -f Makefile  -a  $$i != Test \
			-a $$i != Util -a $$i != UC  \
	                -a $$i != Ana  -a $$i != Ceren \
	                -a $$i != Data  -a $$i != cosmos \
	                -a $$i != Misc  -a $$i != MuUg \
	                -a $$i != Obso  -a $$i != TeX \
	                -a $$i != lib  -a $$i != UserHook \
	                -a $$i != Doc  -a $$i != Betterone \
	                -a $$i != Scrpt -a $$i != NotForExport ] ; then \
		            $(MAKE) $(TARGET) \
			;fi) \
		; fi \
	done
###
$(LIBRARY)(cBrCSamp.o): $(COSMOSINC)/Zglobalc.h \
	$(COSMOSINC)/Zmass.h \
	$(COSMOSINC)/ZbpCnst.h
$(LIBRARY)(cBrLSamp.o): $(COSMOSINC)/ZbpCnst.h \
	$(COSMOSINC)/ZbpTable.h \
	$(COSMOSINC)/Zmass.h
$(LIBRARY)(cBrSSamp.o): $(COSMOSINC)/ZbpCnst.h \
	$(COSMOSINC)/ZbpTable.h \
	$(COSMOSINC)/Zmass.h
$(LIBRARY)(cKnock.o): $(COSMOSINC)/Zptcl.h \
	$(COSMOSINC)/Zunionmap.h \
	$(COSMOSINC)/Zmass.h \
	$(COSMOSINC)/ZdedxAir.h
$(LIBRARY)(cPrCSamp.o): $(COSMOSINC)/ZbpCnst.h \
	$(COSMOSINC)/Zmass.h
$(LIBRARY)(cPrLSamp.o): $(COSMOSINC)/ZbpCnst.h \
	$(COSMOSINC)/ZbpTable.h
$(LIBRARY)(cPrTSamp.o): $(COSMOSINC)/Zmass.h \
	$(COSMOSINC)/ZbpCnst.h
$(LIBRARY)(canihiPath.o): $(COSMOSINC)/Zmass.h
$(LIBRARY)(cbhabhaPath.o): $(COSMOSINC)/Zmass.h
$(LIBRARY)(cbremErgLPM.o): $(COSMOSINC)/Zmass.h
$(LIBRARY)(cbremsPath.o): $(COSMOSINC)/ZbpCnst.h \
	$(COSMOSINC)/ZbpTable.h
$(LIBRARY)(ccomptPath.o): $(COSMOSINC)/Zmass.h
$(LIBRARY)(cdedx2atLE.o): $(COSMOSINC)/ZdedxAir.h
$(LIBRARY)(cdedxElec.o): $(COSMOSINC)/ZdedxAir.h
$(LIBRARY)(cdedxEleci.o): $(COSMOSINC)/ZdedxAir.h
$(LIBRARY)(cdedxInAir.o): $(COSMOSINC)/Zptcl.h \
	$(COSMOSINC)/Zunionmap.h \
	$(COSMOSINC)/Zcode.h \
	$(COSMOSINC)/Zelemagp.h
$(LIBRARY)(cdedxNone.o): $(COSMOSINC)/Zptcl.h \
	$(COSMOSINC)/Zunionmap.h \
	$(COSMOSINC)/Zmass.h \
	$(COSMOSINC)/ZdedxAir.h
$(LIBRARY)(cdedxatHE.o): $(COSMOSINC)/Zptcl.h \
	$(COSMOSINC)/Zunionmap.h \
	$(COSMOSINC)/ZdedxAir.h
$(LIBRARY)(cdedxe.o): $(COSMOSINC)/Zptcl.h \
	$(COSMOSINC)/Zunionmap.h \
	$(COSMOSINC)/ZdedxAir.h \
	$(COSMOSINC)/Zmass.h
$(LIBRARY)(cexcessLen.o): $(COSMOSINC)/Ztrack.h \
	$(COSMOSINC)/Zcondc.h \
	$(COSMOSINC)/Zptcl.h \
	$(COSMOSINC)/Zunionmap.h \
	$(COSMOSINC)/Zcoord.h \
	$(COSMOSINC)/Zpos.h \
	$(COSMOSINC)/Zdirec.h \
	$(COSMOSINC)/Zmagfield.h \
	$(COSMOSINC)/Ztrackv.h \
	$(COSMOSINC)/Zmaxdef.h \
	$(COSMOSINC)/Zelemagp.h
$(LIBRARY)(cmigdFunc.o): $(COSMOSINC)/Zmass.h
$(LIBRARY)(cmollerPath.o): $(COSMOSINC)/Zmass.h
$(LIBRARY)(cpairErgLPM.o): $(COSMOSINC)/Zmass.h
$(LIBRARY)(cpairPath.o): $(COSMOSINC)/ZbpCnst.h
$(LIBRARY)(cphotoEEPath.o): $(COSMOSINC)/ZbpCnst.h \
	$(COSMOSINC)/Zmass.h
$(LIBRARY)(cSampMol.o): $(COSMOSINC)/Ztrackv.h \
	$(COSMOSINC)/Zglobalc.h \
	$(COSMOSINC)/Zmass.h \
	$(COSMOSINC)/Zcode.h

$(LIBRARY):	\
		$(LIBRARY)(cBrCSamp.o) \
		$(LIBRARY)(cBrLSamp.o) \
		$(LIBRARY)(cBrSSamp.o) \
		$(LIBRARY)(cSampMolRedA.o) \
		$(LIBRARY)(cGetXXsec.o) \
               $(LIBRARY)(cKnock.o) \
               $(LIBRARY)(cPrCSamp.o) \
               $(LIBRARY)(cPrLSamp.o) \
               $(LIBRARY)(cPrTSamp.o) \
               $(LIBRARY)(canihiPath.o) \
               $(LIBRARY)(cbhabhaPath.o) \
               $(LIBRARY)(cbremErgLPM.o) \
               $(LIBRARY)(cbremLMPXsec.o) \
               $(LIBRARY)(cbremsPath.o) \
               $(LIBRARY)(ccomptPath.o) \
               $(LIBRARY)(cdedx2atLE.o) \
               $(LIBRARY)(cdedxElec.o) \
               $(LIBRARY)(cdedxEleci.o) \
               $(LIBRARY)(cdedxInAir.o) \
               $(LIBRARY)(cdedxNone.o) \
               $(LIBRARY)(cdedxatHE.o) \
               $(LIBRARY)(cdedxdlt.o) \
               $(LIBRARY)(cdedxe.o) \
               $(LIBRARY)(cexcessLen.o) \
               $(LIBRARY)(cmigdFunc.o) \
               $(LIBRARY)(cmoliere.o) \
               $(LIBRARY)(cmollerPath.o) \
               $(LIBRARY)(cpairErgLPM.o) \
               $(LIBRARY)(cpairLPMXsec.o) \
               $(LIBRARY)(cpairPath.o) \
               $(LIBRARY)(cphotoEEPath.o)  \
               $(LIBRARY)(cMolfunc.o)  \
               $(LIBRARY)(cSampMol.o)
