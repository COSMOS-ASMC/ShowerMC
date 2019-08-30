include  $(LIBLOFT)/site.config


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
		$(LIBLOFTINC)/ZdedxAir.h \
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
		$(LIBLOFTINC)/Ztrack.h \
		$(LIBLOFTINC)/Ztrackp.h \
		$(LIBLOFTINC)/Ztrackv.h \
		$(LIBLOFTINC)/Zunionmap.h \
		$(LIBLOFTINC)/Zxsectionp.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = c2lowerCase.o \
		ckf2cos.o \
		cerrorMsg.o \
		cgetLoginN.o \
		cmanager.o \
		cmkSeed.o \
		copenf.o \
		creadParam.o \
		cskipComment.o \
		cwriteParam.o \
		cwriteSeed.o 


PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = c2lowerCase.f \
		ckf2cos.f \
		cerrorMsg.f \
		cgetLoginN.f \
		cmkSeed.f \
		copenf.f \
		creadParam.f \
		cskipComment.f \
		cwriteParam.f \
		cwriteSeed.f 

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
		cp -r  $(LIBLOFT)/Data  $(LIBLOFT)/Util *
		 $(LIBLOFT)/UserHook  $(LIBLOFT)/Scrpt *
		 $(LIBLOFT)/Doc        $(TEMPLATE)/
		cp $(LIBLOFT)/site.config $(TEMPLATE)/
		@echo "****" Installing header files into $(HEADER)
		cp -r $(LIBLOFT)/Header  $(HEADER)/
		@echo "****" Installation ended.
		@echo "****" Owner and group should be changed by hand

print:
	@$(PRINT) $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"

tags:           $(HDRS) $(SRCS)
	@ctags $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"


Checkdir:
	@( cd $(LIBLOFT)/lib; if [ ! -d  $(ARCH) ]; then  mkdir $(ARCH); fi)



SeeSubDir:
	@for i in  * \
	;do   \
		if [ -d  $$i ] ; then \
			(cd $$i; \
			if [ -f Makefile  -a  $$i != Test \
			-a $$i != Util -a $$i != UC  \
	                -a $$i != Ana  -a $$i != Ceren \
	                -a $$i != Data  -a $$i != Header \
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
$(LIBRARY)(ckf2cos.o): $(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zkfcode.h \
	$(LIBLOFTINC)/Zptcl.h 
$(LIBRARY)(cerrorMsg.o): $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h
$(LIBRARY)(cmanager.o): $(LIBLOFTINC)/Zevhnp.h
$(LIBRARY)(copenf.o): $(LIBLOFTINC)/Zreadonly.h \
	$(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h
$(LIBRARY)(creadParam.o): $(LIBLOFTINC)/ZincForNameL.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zheavyp.h \
	$(LIBLOFTINC)/Zelemagp.h \
	$(LIBLOFTINC)/Zxsectionp.h \
	$(LIBLOFTINC)/Zevhnp.h \
	$(LIBLOFTINC)/Zincidentp.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Znamelist.h \
	$(LIBLOFTINC)/Znamelist1.h \
	$(LIBLOFTINC)/Znamelist2.h
$(LIBRARY)(cwriteParam.o): $(LIBLOFTINC)/ZincForNameL.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zheavyp.h \
	$(LIBLOFTINC)/Zelemagp.h \
	$(LIBLOFTINC)/Zxsectionp.h \
	$(LIBLOFTINC)/Zevhnp.h \
	$(LIBLOFTINC)/Zincidentp.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Znamelist.h \
	$(LIBLOFTINC)/Znamelist1.h \
	$(LIBLOFTINC)/Znamelist2.h
$(LIBRARY)(cwriteSeed.o): $(LIBLOFTINC)/Zmanager.h \
	$(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h

$(LIBRARY):         \
		$(LIBRARY)(c2lowerCase.o) \
		$(LIBRARY)(ckf2cos.o) \
		$(LIBRARY)(cerrorMsg.o) \
		$(LIBRARY)(cgetLoginN.o) \
               $(LIBRARY)(cmkSeed.o) \
               $(LIBRARY)(copenf.o) \
               $(LIBRARY)(creadParam.o) \
               $(LIBRARY)(cskipComment.o) \
               $(LIBRARY)(cwriteParam.o) \
               $(LIBRARY)(cwriteSeed.o)  

