include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcoord.h \
		$(LIBLOFTINC)/Zglobalc.h \
		$(LIBLOFTINC)/Zmagfield.h \
		$(LIBLOFTINC)/Zmass.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/Zunionmap.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = cgUpsilon.o \
		cgetBsin.o \
		ck23.o \
		cmBremE.o \
		cmBremE1.o \
		cmBremE2.o \
		cmBremF1.o \
		cmBremF2.o \
		cmBremI1.o \
		cmBremI2.o \
		cmBremMFP.o \
		cmPairE.o \
		cmPairMFP.o \
		cmPairSpec.o \
		csyncTELoss.o \
		ctXai.o \
		cupsilon.o \
		cxai.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = cgUpsilon.f \
		cgetBsin.f \
		ck23.f \
		cmBremE.f \
		cmBremE1.f \
		cmBremE2.f \
		cmBremF1.f \
		cmBremF2.f \
		cmBremI1.f \
		cmBremI2.f \
		cmBremMFP.f \
		cmPairE.f \
		cmPairMFP.f \
		cmPairSpec.f \
		csyncTELoss.f \
		ctXai.f \
		cupsilon.f \
		cxai.f

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
$(LIBRARY)(cgetBsin.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zmagfield.h
$(LIBRARY)(cmBremMFP.o): $(LIBLOFTINC)/Zglobalc.h
$(LIBRARY)(cmPairMFP.o): $(LIBLOFTINC)/Zglobalc.h
$(LIBRARY)(csyncTELoss.o): $(LIBLOFTINC)/Zglobalc.h
$(LIBRARY)(cupsilon.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmagfield.h
$(LIBRARY)(cxai.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmass.h \
	$(LIBLOFTINC)/Zmagfield.h

$(LIBRARY):         \
               $(LIBRARY)(cgUpsilon.o) \
               $(LIBRARY)(cgetBsin.o) \
               $(LIBRARY)(ck23.o) \
               $(LIBRARY)(cmBremE.o) \
               $(LIBRARY)(cmBremE1.o) \
               $(LIBRARY)(cmBremE2.o) \
               $(LIBRARY)(cmBremF1.o) \
               $(LIBRARY)(cmBremF2.o) \
               $(LIBRARY)(cmBremI1.o) \
               $(LIBRARY)(cmBremI2.o) \
               $(LIBRARY)(cmBremMFP.o) \
               $(LIBRARY)(cmPairE.o) \
               $(LIBRARY)(cmPairMFP.o) \
               $(LIBRARY)(cmPairSpec.o) \
               $(LIBRARY)(csyncTELoss.o) \
               $(LIBRARY)(ctXai.o) \
               $(LIBRARY)(cupsilon.o) \
               $(LIBRARY)(cxai.o) 
