include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zmanagerp.h \
		$(LIBLOFTINC)/Zmaxdef.h \
		$(LIBLOFTINC)/Zprimary.h \
		$(LIBLOFTINC)/Zprimaryc.h \
		$(LIBLOFTINC)/Zprimaryv.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/ZrigCut.h \
		$(LIBLOFTINC)/Zunionmap.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = ciniSPrim.o \
		cmkPrimSTbl.o \
		cprintPrim.o \
		cprocPrimDt.o \
		crdPrimData.o \
		crigCut.o \
		csampPrimary.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = ciniSPrim.f \
		cmkPrimSTbl.f \
		cprintPrim.f \
		cprocPrimDt.f \
		crdPrimData.f \
		crigCut.f \
		csampPrimary.f

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
$(LIBRARY)(ciniSPrim.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zprimary.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zprimaryc.h \
	$(LIBLOFTINC)/Zprimaryv.h \
	$(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cmkPrimSTbl.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zprimary.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zprimaryc.h
$(LIBRARY)(cprintPrim.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zprimary.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zprimaryv.h
$(LIBRARY)(cprocPrimDt.o): $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zprimary.h \
	$(LIBLOFTINC)/Zprimaryc.h
$(LIBRARY)(crdPrimData.o): $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zprimary.h \
	$(LIBLOFTINC)/Zprimaryc.h \
	$(LIBLOFTINC)/Zprimaryv.h
$(LIBRARY)(crigCut.o): $(LIBLOFTINC)/ZrigCut.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zmanagerp.h
$(LIBRARY)(csampPrimary.o): $(LIBLOFTINC)/Zmanagerp.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zprimary.h \
	$(LIBLOFTINC)/Zprimaryv.h \
	$(LIBLOFTINC)/Zcode.h

$(LIBRARY):         \
               $(LIBRARY)(ciniSPrim.o) \
               $(LIBRARY)(cmkPrimSTbl.o) \
               $(LIBRARY)(cprintPrim.o) \
               $(LIBRARY)(cprocPrimDt.o) \
               $(LIBRARY)(crdPrimData.o) \
               $(LIBRARY)(crigCut.o) \
               $(LIBRARY)(csampPrimary.o) 
