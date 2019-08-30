include  $(COSMOSTOP)/site.config


EXTHDRS	      = $(LIBLOFTINC)/BlockData/cblkElemag.h \
		$(LIBLOFTINC)/BlockData/cblkEvhnp.h \
		$(LIBLOFTINC)/BlockData/cblkGene.h \
		$(LIBLOFTINC)/BlockData/cblkHeavy.h \
		$(LIBLOFTINC)/BlockData/cblkIncident.h \
		$(LIBLOFTINC)/BlockData/cblkManager.h \
		$(LIBLOFTINC)/BlockData/cblkMuInt.h \
		$(LIBLOFTINC)/BlockData/cblkObs.h \
		$(LIBLOFTINC)/BlockData/cblkSpecial.h \
		$(LIBLOFTINC)/BlockData/cblkTracking.h \
		$(LIBLOFTINC)/BlockData/cblkXsec.h \
		$(LIBLOFTINC)/BlockData/cblkdedx.h \
		$(LIBLOFTINC)/Zcmuint.h \
		$(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zcondc.h \
		$(LIBLOFTINC)/Zcoord.h \
		$(LIBLOFTINC)/ZcosmosBD.h \
		$(COSMOSINC)/ZdedxAir.h \
		$(LIBLOFTINC)/Zdirec.h \
		$(LIBLOFTINC)/Zelemagp.h \
		$(LIBLOFTINC)/Zevhnp.h \
		$(LIBLOFTINC)/Zheavyp.h \
		$(LIBLOFTINC)/Zheavyv.h \
		$(LIBLOFTINC)/ZincForNameL.h \
		$(LIBLOFTINC)/Zincidentp.h \
		$(LIBLOFTINC)/Zmagfield.h \
		$(LIBLOFTINC)/Zmanager.h \
		$(LIBLOFTINC)/Zmanagerp.h \
		$(LIBLOFTINC)/Zmaxdef.h \
		$(LIBLOFTINC)/Znamelist.h \
		$(LIBLOFTINC)/Znamelist1.h \
		$(LIBLOFTINC)/Znamelist2.h \
		$(LIBLOFTINC)/Zobs.h \
		$(LIBLOFTINC)/Zobsp.h \
		$(LIBLOFTINC)/Zobsv.h \
		$(LIBLOFTINC)/Zpos.h \
		$(LIBLOFTINC)/Zprimary.h \
		$(LIBLOFTINC)/Zprimaryv.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/Zreadonly.h \
		$(LIBLOFTINC)/Zspecial.h \
		$(COSMOSINC)/Ztrack.h \
		$(LIBLOFTINC)/Ztrackp.h \
		$(COSMOSINC)/Ztrackv.h \
		$(LIBLOFTINC)/Zunionmap.h \
		$(LIBLOFTINC)/Zxsectionp.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = cbeginRun.o \
		cendRun.o \
		ceventLoop.o \
		cmanager.o \
		cqEventNo.o \
		cqUHookr.o \
		cqverstion.o


PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = cbeginRun.f \
		cendRun.f \
		ceventLoop.f \
		cmanager.f \
		cqEventNo.f \
		cqUHookr.f \
		cqverstion.f


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
$(LIBRARY)(cbeginRun.o): $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zmanager.h \
	$(LIBLOFTINC)/Zelemagp.h \
	$(LIBLOFTINC)/Zevhnp.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(COSMOSINC)/Ztrackv.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(LIBLOFTINC)/Zincidentp.h \
	$(LIBLOFTINC)/Zprimary.h \
	$(LIBLOFTINC)/Zprimaryv.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/Zheavyv.h
$(LIBRARY)(cendRun.o): $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zprimary.h \
	$(LIBLOFTINC)/Zprimaryv.h

$(LIBRARY)(ceventLoop.o): $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zmanager.h \
	$(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zprimary.h \
	$(LIBLOFTINC)/Zprimaryv.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zobsv.h
$(LIBRARY)(cmanager.o): $(LIBLOFTINC)/Zevhnp.h
$(LIBRARY)(cqEventNo.o): $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h
$(LIBRARY)(cqUHookr.o): $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h

$(LIBRARY):         \
               $(LIBRARY)(cbeginRun.o) \
               $(LIBRARY)(cendRun.o) \
               $(LIBRARY)(ceventLoop.o) \
               $(LIBRARY)(cmanager.o) \
               $(LIBRARY)(cqEventNo.o) \
               $(LIBRARY)(cqUHookr.o) \
               $(LIBRARY)(cqverstion.o)
