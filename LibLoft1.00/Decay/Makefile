include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/Zunionmap.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = cbomegaDcy.o \
		cdDecay.o \
		cetaDecay.o \
		cgzaiDecay.o \
		ckaonDecay.o \
		clambdaDcy.o \
		clambdacDcy.o \
		cmuNeuDcy.o \
		cpCos2pxyz.o \
		cpi0Decay.o \
		cpiMuDecay.o \
		csigmaDecay.o \
		cvmDecay.o \
		ctauNeuDcy.o \
		cetapDecay.o \
		cDsDecay.o \
		cDeltaDecay.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = cbomegaDcy.f \
		cdDecay.f \
		cetaDecay.f \
		cgzaiDecay.f \
		ckaonDecay.f \
		clambdaDcy.f \
		clambdacDcy.f \
		cmuNeuDcy.f \
		cpCos2pxyz.f \
		cpi0Decay.f \
		cpiMuDecay.f \
		csigmaDecay.f \
		cvmDecay.f \
		ctauNeuDcy.f  \
		cetapDecay.f \
		cDsDecay.f \
		cDeltaDecay.f

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
		cp -r $(LIBLOFT)/cosmos  $(HEADER)/
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
$(LIBRARY)(cbomegaDcy.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cdDecay.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cetaDecay.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cgzaiDecay.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(ckaonDecay.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(clambdaDcy.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(clambdacDcy.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cmuNeuDcy.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cpCos2pxyz.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cpi0Decay.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cpiMuDecay.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(csigmaDecay.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cvmDecay.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h

$(LIBRARY):         \
               $(LIBRARY)(cbomegaDcy.o) \
               $(LIBRARY)(cdDecay.o) \
               $(LIBRARY)(cetaDecay.o) \
               $(LIBRARY)(cgzaiDecay.o) \
               $(LIBRARY)(ckaonDecay.o) \
               $(LIBRARY)(clambdaDcy.o) \
               $(LIBRARY)(clambdacDcy.o) \
               $(LIBRARY)(cmuNeuDcy.o) \
               $(LIBRARY)(cpCos2pxyz.o) \
               $(LIBRARY)(cpi0Decay.o) \
               $(LIBRARY)(cpiMuDecay.o) \
               $(LIBRARY)(csigmaDecay.o) \
               $(LIBRARY)(cvmDecay.o) \
               $(LIBRARY)(ctauNeuDcy.o) \
               $(LIBRARY)(cetapDecay.o) \
               $(LIBRARY)(cDsDecay.o) \
               $(LIBRARY)(cDeltaDecay.o)
