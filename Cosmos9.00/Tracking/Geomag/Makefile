include  $(COSMOSTOP)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zcondc.h \
		$(LIBLOFTINC)/Zcoord.h \
		$(LIBLOFTINC)/Zdirec.h \
		$(LIBLOFTINC)/Zglobalc.h \
		$(LIBLOFTINC)/Zmagfield.h \
		$(LIBLOFTINC)/Zmagfieldtp.h \
		$(LIBLOFTINC)/Zmaxdef.h \
		$(LIBLOFTINC)/Zobs.h \
		$(LIBLOFTINC)/Zobsp.h \
		$(LIBLOFTINC)/Zpos.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(COSMOSINC)/Ztrack.h \
		$(LIBLOFTINC)/Zunionmap.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = cbDefByRK.o \
		cgeomag.o \
		cmagneticDef.o \
		cprintMagF.o \
		csetMagField.o \
		ctransMagTo.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = cbDefByRK.f \
		cgeomag.f \
		cmagneticDef.f \
		cprintMagF.f \
		csetMagField.f \
		ctransMagTo.f

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
$(LIBRARY)(cbDefByRK.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zobsp.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cgeomag.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmagfield.h 
$(LIBRARY)(cmagneticDef.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h
$(LIBRARY)(cprintMagF.o): $(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(csetMagField.o): $(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmagfieldtp.h
$(LIBRARY)(ctransMagTo.o): $(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zglobalc.h

$(LIBRARY):         \
               $(LIBRARY)(cbDefByRK.o) \
               $(LIBRARY)(cgeomag.o) \
               $(LIBRARY)(cmagneticDef.o) \
               $(LIBRARY)(cprintMagF.o) \
               $(LIBRARY)(csetMagField.o) \
               $(LIBRARY)(ctransMagTo.o) 
