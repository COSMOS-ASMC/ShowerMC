include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcmuint.h \
		$(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zcondc.h \
		$(LIBLOFTINC)/Zcoord.h \
		$(LIBLOFTINC)/Zdirec.h \
		$(LIBLOFTINC)/Zmagfield.h \
		$(LIBLOFTINC)/Zmass.h \
		$(LIBLOFTINC)/Zmaxdef.h \
		$(LIBLOFTINC)/Zpos.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/Ztrack.h \
		$(LIBLOFTINC)/Ztrackp.h \
		$(LIBLOFTINC)/Ztrackv.h \
		$(LIBLOFTINC)/Zunionmap.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = cmuBrsmp.o \
		cmuInte.o \
		cmuNsmp.o \
		cmuPrsmp.o \
		cmudEdx.o \
		cmuvmax2.o \
		crdmuSmpTbl.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = cmuBrsmp.f \
		cmuInte.f \
		cmuNsmp.f \
		cmuPrsmp.f \
		cmudEdx.f \
		cmuvmax2.f \
		crdmuSmpTbl.f

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
$(LIBRARY)(cmuBrsmp.o): $(LIBLOFTINC)/Zcmuint.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(cmuInte.o): $(LIBLOFTINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Ztrackv.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Ztrackp.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cmuNsmp.o): $(LIBLOFTINC)/Zcmuint.h
$(LIBRARY)(cmuPrsmp.o): $(LIBLOFTINC)/Zcmuint.h
$(LIBRARY)(cmudEdx.o): $(LIBLOFTINC)/Zcmuint.h
$(LIBRARY)(cmuvmax2.o): $(LIBLOFTINC)/Zcmuint.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(crdmuSmpTbl.o): $(LIBLOFTINC)/Zcmuint.h

$(LIBRARY):         \
               $(LIBRARY)(cmuBrsmp.o) \
               $(LIBRARY)(cmuInte.o) \
               $(LIBRARY)(cmuNsmp.o) \
               $(LIBRARY)(cmuPrsmp.o) \
               $(LIBRARY)(cmudEdx.o) \
               $(LIBRARY)(cmuvmax2.o) \
               $(LIBRARY)(crdmuSmpTbl.o) 
