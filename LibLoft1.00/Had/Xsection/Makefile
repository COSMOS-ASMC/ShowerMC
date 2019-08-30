include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/Zunionmap.h \
		$(LIBLOFTINC)/Zxsectionp.h \
		$(LIBLOFTINC)/Ztrackp.h  \
		$(LIBLOFTINC)/Ztrack.h  \
		$(LIBLOFTINC)/Ztrackv.h 

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = cpdgXs.o \
		cAAXsec2.o \
		cgpXsec.o \
		cgpxs1.o \
		cppXsec.o \
	        cprotonAXsec.o \
		ctotx.o \
		cxp2xAXsec.o \
		cxAbyxpXsec.o \
		ckmnTotXs.o \
		ckmpTotXs.o \
		ckpnTotXs.o \
		ckppTotXs.o \
		cnbarpTotXs.o \
		cnpTotXs.o \
		cpbarpTotXs.o \
		cpimnTotXs.o \
		cpimpTotXs.o \
		cpipnTotXs.o \
		cpippTotXs.o \
		cppTotXs.o  \
		cppXsec.o \
	        cprotonAXsec.o \
		cinelx.o \
		cpdgXs200GeV.o \
		cGetXsec.o \
		cxsSpecial.o


PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = cpdgXs.f \
		cAAXsec2.f \
		cgpXsec.f \
		cgpxs1.f \
		ctotx.f \
		cxp2xAXsec.f \
		cxAbyxpXsec.f \
		ckmnTotXs.f \
		ckmpTotXs.f \
		ckpnTotXs.f \
		ckppTotXs.f \
		cnbarpTotXs.f \
		cnpTotXs.f \
		cpbarpTotXs.f \
		cpimnTotXs.f \
		cpimpTotXs.f \
		cpipnTotXs.f \
		cpippTotXs.f \
		cppTotXs.f  \
		cppXsec.f \
	        cprotonAXsec.f \
		cinelx.f \
		cpdgXs200GeV.f90  \
		cGetXsec.f \
		cxsSpecial.f \

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
$(LIBRARY)(cpdgXs.o): $(LIBLOFTINC)/Zmass.h 
$(LIBRARY)(cAAXsec2.o): $(LIBLOFTINC)/Zxsectionp.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zptcl.h
$(LIBRARY)(cppXsec.o): $(LIBLOFTINC)/Zxsectionp.h

$(LIBRARY)(ctotx.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h

$(LIBRARY)(cxsSpecial.o): $(LIBLOFTINC)/Zcode.h \
        $(LIBLOFTINC)/Zptcl.h \
        $(LIBLOFTINC)/Zevhnp.h
$(LIBRARY):         \
		$(LIBRARY)(cpdgXs.o) \
		$(LIBRARY)(cAAXsec2.o) \
		$(LIBRARY)(cgpXsec.o) \
		$(LIBRARY)(cgpxs1.o) \
		$(LIBRARY)(ctotx.o) \
		$(LIBRARY)(cxp2xAXsec.o) \
		$(LIBRARY)(cxAbyxpXsec.o) \
		$(LIBRARY)(ckmnTotXs.o) \
		$(LIBRARY)(ckmpTotXs.o) \
		$(LIBRARY)(ckpnTotXs.o) \
		$(LIBRARY)(ckppTotXs.o) \
		$(LIBRARY)(cnbarpTotXs.o) \
		$(LIBRARY)(cnpTotXs.o) \
		$(LIBRARY)(cpbarpTotXs.o) \
		$(LIBRARY)(cpimnTotXs.o) \
		$(LIBRARY)(cpimpTotXs.o) \
		$(LIBRARY)(cpipnTotXs.o) \
		$(LIBRARY)(cpippTotXs.o) \
		$(LIBRARY)(cppTotXs.o)  \
		$(LIBRARY)(cppXsec.o) \
		$(LIBRARY)(cprotonAXsec.o) \
		$(LIBRARY)(cinelx.o) \
		$(LIBRARY)(cpdgXs200GeV.o) \
		$(LIBRARY)(cGetXsec.o) \
		$(LIBRARY)(cxsSpecial.o)
