include  $(LIBLOFT)/site.config



EXTHDRS	      = $(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zcoord.h \
		$(LIBLOFTINC)/Zdirec.h \
		$(LIBLOFTINC)/Zelemagp.h \
		$(LIBLOFTINC)/Zevhnp.h \
		$(LIBLOFTINC)/Zglobalc.h \
		$(LIBLOFTINC)/Zptcl.h


HDRS	      =


MAKEFILE      = Makefile

OBJS	      = cepControl.o \
		epprocAlias.o \
		epSrim.o \
		cHowMCS.o \
		csampAF.o  \
		epSmpBremAng2BN.o \
		epSmpBremAngRelF.o \
		epSmpBremAngTsai.o \
		kggamma.o  \
		cNucMass.o \
		cTargetFrag.o \
		cCodeConf.o


PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = cepControl.f \
		epprocAlias.f \
		epSrim.f90 \
		cHowMCS.f90 \
		csampAF.f90  \
		epSmpBremAng2BN.f90 \
		epSmpBremAngRelF.f90 \
		epSmpBremAngTsai.f90 \
		kggamma.f90  \
		cNucMass.f90 \
		cTargetFrag.f90 \
		cCodeConv.f

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
$(LIBRARY)(cepControl.o): $(LIBLOFTINC)/Zptcl.h 


$(LIBRARY):         \
               $(LIBRARY)(cepControl.o) \
	       $(LIBRARY)(epprocAlias.o) \
	       $(LIBRARY)(epSrim.o) \
	       $(LIBRARY)(cHowMCS.o) \
	       $(LIBRARY)(csampAF.o) \
	       $(LIBRARY)(epSmpBremAng2BN.o) \
	       $(LIBRARY)(epSmpBremAngRelF.o) \
	       $(LIBRARY)(epSmpBremAngTsai.o) \
	       $(LIBRARY)(kggamma.o) \
	       $(LIBRARY)(cNucMass.o) \
		$(LIBRARY)(cTargetFrag.o) \
		$(LIBRARY)(cCodeConv.o) 
