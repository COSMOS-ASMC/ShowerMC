include  $(COSMOSTOP)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcondc.h \
		$(LIBLOFTINC)/Zcoord.h \
		$(LIBLOFTINC)/Zdirec.h \
		$(LIBLOFTINC)/Zelemagp.h \
		$(LIBLOFTINC)/Zmagfield.h \
		$(LIBLOFTINC)/Zmaxdef.h \
		$(LIBLOFTINC)/Zobs.h \
		$(LIBLOFTINC)/Zobsv.h \
		$(LIBLOFTINC)/Zpos.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(COSMOSINC)/Ztrack.h \
		$(LIBLOFTINC)/Zunionmap.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = cNeByApproxB.o \
		cavedEdx.o \
		cgetMoliereU.o \
		cgetNmu.o \
		clatDist.o \
		cobAS.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = cNeByApproxB.f \
		cavedEdx.f \
		cgetMoliereU.f \
		cgetNmu.f \
		clatDist.f \
		cobAS.f

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
$(LIBRARY)(cgetMoliereU.o): $(LIBLOFTINC)/Zelemagp.h
$(LIBRARY)(cgetNmu.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/Zelemagp.h
$(LIBRARY)(cobAS.o): $(COSMOSINC)/Ztrack.h \
	$(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcoord.h \
	$(LIBLOFTINC)/Zpos.h \
	$(LIBLOFTINC)/Zdirec.h \
	$(LIBLOFTINC)/Zmagfield.h \
	$(LIBLOFTINC)/Zobs.h \
	$(LIBLOFTINC)/Zmaxdef.h \
	$(LIBLOFTINC)/Zobsv.h \
	$(LIBLOFTINC)/Zelemagp.h

$(LIBRARY):         \
               $(LIBRARY)(cNeByApproxB.o) \
               $(LIBRARY)(cavedEdx.o) \
               $(LIBRARY)(cgetMoliereU.o) \
               $(LIBRARY)(cgetNmu.o) \
               $(LIBRARY)(clatDist.o) \
               $(LIBRARY)(cobAS.o) 
