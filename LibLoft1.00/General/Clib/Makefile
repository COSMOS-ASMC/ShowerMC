include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zheavyp.h \
		$(LIBLOFTINC)/Zmass.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/ZsaveStruc.h \
		$(LIBLOFTINC)/Zunionmap.h 


HDRS	      =


MAKEFILE      = Makefile

OBJS	      = c2bdcy.o \
		cadjm.o \
		cbst0.o \
		cbst1.o \
		cgeqm.o \
		cgetRotMat4.o \
		cgetRotMat3.o \
		cgetcm.o \
		cgetBeta.o \
		cgetctau.o \
		cgetlf.o \
		cibst1.o \
		cloreb.o \
		clorep.o \
		clorez.o \
		cmkptc.o \
		cnbdcy.o \
		cpm2e.o \
		cpxyzp.o \
		crot3mom.o \
		crot3vec.o \
		ctestOnShell.o \
		cscalerProd.o \
		cvecProd.o \
		ctransVectZ.o \
		ce2pp.o \
		cyeta.o 



PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = c2bdcy.f \
		cadjm.f \
		cbst0.f \
		cbst1.f \
		cgeqm.f \
		cgetRotMat4.f \
		cgetRotMat3.f \
		cgetcm.f \
		cgetBeta.f \
		cgetctau.f \
		cgetlf.f \
		cibst1.f \
		cloreb.f \
		clorep.f \
		clorez.f90 \
		cmkptc.f \
		cnbdcy.f \
		cpm2e.f \
		cpxyzp.f \
		crot3mom.f \
		crot3vec.f \
		ctestOnShell.f \
		cscalerProd.f \
		cvecProd.f \
		ctransVectZ.f \
		ce2pp.f \
		cyeta.f 



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
			-a $$i != Obso \
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
$(LIBRARY)(c2bdcy.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cadjm.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cbst0.o): $(LIBLOFTINC)/ZsaveStruc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cbst1.o): $(LIBLOFTINC)/ZsaveStruc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cgeqm.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cgetRotMat4.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cgetcm.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cgetlf.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cibst1.o): $(LIBLOFTINC)/ZsaveStruc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cloreb.o): $(LIBLOFTINC)/ZsaveStruc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(clorep.o): $(LIBLOFTINC)/ZsaveStruc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(clorez.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cmkptc.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zheavyp.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(cnbdcy.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cpm2e.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cpxyzp.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(crot3mom.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(crot3vec.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(ctestOnShell.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cyeta.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h


$(LIBRARY):         \
               $(LIBRARY)(c2bdcy.o) \
               $(LIBRARY)(cadjm.o) \
               $(LIBRARY)(cbst0.o) \
               $(LIBRARY)(cbst1.o) \
               $(LIBRARY)(cgeqm.o) \
               $(LIBRARY)(cgetRotMat4.o) \
               $(LIBRARY)(cgetRotMat3.o) \
               $(LIBRARY)(cgetcm.o) \
               $(LIBRARY)(cgetBeta.o) \
               $(LIBRARY)(cgetctau.o) \
               $(LIBRARY)(cgetlf.o) \
               $(LIBRARY)(cibst1.o) \
               $(LIBRARY)(cloreb.o) \
               $(LIBRARY)(clorep.o) \
               $(LIBRARY)(clorez.o) \
               $(LIBRARY)(cmkptc.o) \
               $(LIBRARY)(cnbdcy.o) \
               $(LIBRARY)(cpm2e.o) \
               $(LIBRARY)(cpxyzp.o) \
               $(LIBRARY)(crot3mom.o) \
               $(LIBRARY)(crot3vec.o) \
               $(LIBRARY)(ctestOnShell.o) \
	       $(LIBRARY)(cscalerProd.o) \
	       $(LIBRARY)(cvecProd.o) \
		$(LIBRARY)(ctransVectZ.o) \
		$(LIBRARY)(ce2pp.o) \
               $(LIBRARY)(cyeta.o) 
