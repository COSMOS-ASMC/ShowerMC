include  $(COSMOSTOP)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcmuint.h \
		$(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zcondc.h \
		$(LIBLOFTINC)/Zcoord.h \
		$(COSMOSINC)/Zcoordtype.h \
		$(LIBLOFTINC)/Zdirec.h \
		$(LIBLOFTINC)/Zelemagp.h \
		$(LIBLOFTINC)/Zevhnp.h \
		$(LIBLOFTINC)/Zevhnv.h \
		$(LIBLOFTINC)/Zglobalc.h \
		$(LIBLOFTINC)/Zheavyc.h \
		$(LIBLOFTINC)/Zheavyp.h \
		$(LIBLOFTINC)/Zincidentp.h \
		$(LIBLOFTINC)/Zincidentv.h \
		$(LIBLOFTINC)/Zlife.h \
		$(LIBLOFTINC)/Zmagfield.h \
		$(LIBLOFTINC)/Zmanagerp.h \
		$(LIBLOFTINC)/Zmass.h \
		$(LIBLOFTINC)/Zmaxdef.h \
		$(LIBLOFTINC)/Zobs.h \
		$(LIBLOFTINC)/Zobsp.h \
		$(LIBLOFTINC)/Zobsv.h \
		$(LIBLOFTINC)/Zpos.h \
		$(LIBLOFTINC)/Zprimary.h \
		$(LIBLOFTINC)/Zprimaryv.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/ZrigCut.h \
		$(LIBLOFTINC)/ZsaveStruc.h \
		$(COSMOSINC)/Zstackv.h \
		$(LIBLOFTINC)/ZsubstRec.h \
		$(COSMOSINC)/Ztrack.h \
		$(LIBLOFTINC)/Ztrackp.h \
		$(COSMOSINC)/Ztrackv.h \
		$(LIBLOFTINC)/Zunionmap.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      =  \
		c3DV2DDCos.o \
		c4to3mom.o \
		cdefByMagAndE.o \
		ccompPathEnd.o \
		ce2p.o \
		ceCent2sph.o \
		cfinTracking.o \
		cgetMFP.o \
		cgetZenith.o \
		cifDead.o \
		cifXObsSite.o \
		ciniSmpIntL.o \
		ciniTracking.o \
		cinitObs.o \
		cinteraction.o \
		cmaxMovLen.o \
		cmkIncident.o \
		cobservation.o \
		cprObsSite.o \
		cprintCoord.o \
		cprintObs.o \
		cputCerenkov.o \
		cputTrInfo.o \
		cresetDirec.o \
		cresetMom.o \
		csPrimAng.o \
		csetCoord.o \
		csetPos.o \
		csph2eCent.o \
		cstack.o \
		ctracking.o \
		ctransCoord2.o \
		ctransVectZx.o \
		cvec2dcos.o \
		cxyz2det.o \
		cxyz2prim.o \
		cdet2xyzFE.o \
		cfixProc.o  \
		cpreSampEIntL.o \
		cpreSampGIntL.o \
		cpreSampNEPIntL.o \
		cexcessLen.o  \
		cEfield0.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = \
		c3DV2DDCos.f \
		c4to3mom.f \
		cdefByMagAndE.f \
		ccompPathEnd.f \
		ce2p.f \
		ceCent2sph.f \
		cfinTracking.f \
		cgetMFP.f \
		cgetZenith.f \
		cifDead.f \
		cifXObsSite.f \
		ciniSmpIntL.f \
		ciniTracking.f \
		cinitObs.f \
		cinteraction.f \
		cmaxMovLen.f \
		cmkIncident.f \
		cobservation.f \
		cprObsSite.f \
		cprintCoord.f \
		cprintObs.f \
		cputCerenkov.f \
		cputTrInfo.f \
		cresetDirec.f \
		cresetMom.f \
		csPrimAng.f \
		csetCoord.f \
		csetPos.f \
		csph2eCent.f \
		cstack.f \
		ctracking.f \
		ctransCoord2.f \
		ctransVectZx.f \
		cvec2dcos.f \
		cxyz2det.f \
		cxyz2prim.f \
		cdet2xyzFE.f \
		cfixProc.f \
		cpreSampEIntL.f \
		cpreSampGIntL.f \
		cpreSampNEPIntL.f \
		cexcessLen.f \
		cEfield0.f


SYSHDRS	      =

all:
	make	Checkdir $(LIBRARY) SeeSubDir "TARGET=$@"

clean:
	@rm -f $(OBJS) core  *~ a.out temp*.o 
	make SeeSubDir "TARGET=$@"

veryclean:
	@rm -f $(OBJS) core  *~ a.out temp*.o temp*.f*
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
	@for i in  * \
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
$(LIBRARY)(c3DV2DDCos.o): $(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(c4to3mom.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h
$(LIBRARY)(cdefByMagAndE.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h 


$(LIBRARY)(ccompPathEnd.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zelemagp.h \
	$(LIBLOFTINC)/ZsubstRec.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h
$(LIBRARY)(ce2p.o):  $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h
$(LIBRARY)(ceCent2sph.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cfinTracking.o):   $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h

$(LIBRARY)(cgetZenith.o): $(LIBLOFTINC)/ZsubstRec.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h
$(LIBRARY)(cifDead.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zcode.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h
$(LIBRARY)(cifXObsSite.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zincidentv.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(ciniSmpIntL.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h
$(LIBRARY)(ciniTracking.o): $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/Zincidentv.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zheavyc.h \
	$(LIBLOFTINC)/Zheavyp.h \
	$(LIBLOFTINC)/Zelemagp.h
$(LIBRARY)(cinitObs.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/Zincidentp.h
$(LIBRARY)(cinteraction.o): $(LIBLOFTINC)/ZsubstRec.h \
	$(LIBLOFTINC)/Zcode.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zincidentv.h
$(LIBRARY)(cmaxMovLen.o): $(LIBLOFTINC)/Zcode.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zelemagp.h 

$(LIBRARY)(cmkIncident.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/Zincidentp.h \
	$(LIBLOFTINC)/Zincidentv.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/ZrigCut.h
$(LIBRARY)(cobservation.o): $(LIBLOFTINC)/Zglobalc.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h
$(LIBRARY)(cprObsSite.o): $(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h
$(LIBRARY)(cprintCoord.o): $(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(COSMOSINC)/Zcoordtype.h
$(LIBRARY)(cprintObs.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h
$(LIBRARY)(cputCerenkov.o): $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zmass.h \
	$(LIBLOFTINC)/Zprimary.h \
	$(LIBLOFTINC)/Zprimaryv.h \
	$(LIBLOFTINC)/Zheavyp.h \
	$(LIBLOFTINC)/Zincidentv.h
$(LIBRARY)(cputTrInfo.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsv.h
$(LIBRARY)(cresetDirec.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h
$(LIBRARY)(cresetMom.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h
$(LIBRARY)(csPrimAng.o): $(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zincidentp.h \
	$(LIBLOFTINC)/Zglobalc.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zincidentv.h
$(LIBRARY)(csetCoord.o): $(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(COSMOSINC)/Zcoordtype.h
$(LIBRARY)(csetPos.o): $(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zpos.h 

$(LIBRARY)(csph2eCent.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cstack.o): $(LIBLOFTINC)/Zmaxdef.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(COSMOSINC)/Zstackv.h
$(LIBRARY)(ctracking.o): $(LIBLOFTINC)/Zcode.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/Zincidentv.h \
	$(LIBLOFTINC)/Zmanagerp.h \
	$(COSMOSINC)/Zstackv.h \
	$(LIBLOFTINC)/Zglobalc.h 

$(LIBRARY)(ctransCoord2.o): $(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(ctransVectZx.o): $(LIBLOFTINC)/ZsaveStruc.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cvec2dcos.o): $(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cxyz2det.o): $(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zobsv.h
$(LIBRARY)(cxyz2prim.o): $(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zobsv.h



$(LIBRARY):         \
               $(LIBRARY)(c3DV2DDCos.o) \
               $(LIBRARY)(c4to3mom.o) \
               $(LIBRARY)(cdefByMagAndE.o) \
               $(LIBRARY)(ccompPathEnd.o) \
               $(LIBRARY)(ce2p.o) \
               $(LIBRARY)(ceCent2sph.o) \
               $(LIBRARY)(cfinTracking.o) \
               $(LIBRARY)(cgetMFP.o) \
               $(LIBRARY)(cgetZenith.o) \
               $(LIBRARY)(cifDead.o) \
               $(LIBRARY)(cifXObsSite.o) \
               $(LIBRARY)(ciniSmpIntL.o) \
               $(LIBRARY)(ciniTracking.o) \
               $(LIBRARY)(cinitObs.o) \
               $(LIBRARY)(cinteraction.o) \
               $(LIBRARY)(cmaxMovLen.o) \
               $(LIBRARY)(cmkIncident.o) \
               $(LIBRARY)(cobservation.o) \
               $(LIBRARY)(cprObsSite.o) \
               $(LIBRARY)(cprintCoord.o) \
               $(LIBRARY)(cprintObs.o) \
               $(LIBRARY)(cputCerenkov.o) \
               $(LIBRARY)(cputTrInfo.o) \
               $(LIBRARY)(cresetDirec.o) \
               $(LIBRARY)(cresetMom.o) \
               $(LIBRARY)(csPrimAng.o) \
               $(LIBRARY)(csetCoord.o) \
               $(LIBRARY)(csetPos.o) \
               $(LIBRARY)(csph2eCent.o) \
               $(LIBRARY)(cstack.o) \
               $(LIBRARY)(ctracking.o) \
               $(LIBRARY)(ctransCoord2.o) \
               $(LIBRARY)(ctransVectZx.o) \
               $(LIBRARY)(cvec2dcos.o) \
               $(LIBRARY)(cxyz2det.o) \
               $(LIBRARY)(cxyz2prim.o) \
	       $(LIBRARY)(cdet2xyzFE.o) \
	       $(LIBRARY)(cfixProc.o) \
	       $(LIBRARY)(cpreSampEIntL.o) \
	       $(LIBRARY)(cpreSampGIntL.o) \
	       $(LIBRARY)(cpreSampNEPIntL.o) \
	       $(LIBRARY)(cexcessLen.o) \
	       $(LIBRARY)(cEfield0.o) 


