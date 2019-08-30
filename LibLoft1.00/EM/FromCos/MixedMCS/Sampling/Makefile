include  $(LIBLOFT)/site.config


EXTHDRS       = $(LIBLOFTINC)/Zmanager.h \
                $(LIBLOFTINC)/Zmanagerp.h \
                $(LIBLOFTINC)/Zcode.h \
                $(LIBLOFTINC)/Zptcl.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = cMCS0.o \
		cMCSTypeDef.o \
		cmkMCSconst.o \
		cgetLambda.o \
		ciniMixedMS.o \
		cgetMuc.o \
		csampCS.o \
		csampSoftMCS.o \
		cMCSconstForMedia.o \
		cfixMCSmodel.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = cMCS0.f90  \
		cMCSTypeDef.f90 \
		cgetLambda.f90 \
		cmkMCSconst.f90 \
		ciniMixedMS.f90 \
		cgetMuc.f90  \
		csampCS.f90  \
		csampSoftMCS.f90 \
		cMCSconstForMedia.f \
		cfixMCSmodel.f
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
	@for i  in  * \
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
$(LIBRARY)(cMCSconstForMedia.o): $(LIBLOFTINC)/Zmanager.h \
		$(LIBLOFTINC)/Zmanagerp.h 

$(LIBRARY):	\
		$(LIBRARY)(cMCS0.o) \
		$(LIBRARY)(cMCSTypeDef.o) \
		$(LIBRARY)(cmkMCSconst.o) \
		$(LIBRARY)(cgetLambda.o) \
		$(LIBRARY)(ciniMixedMS.o) \
		$(LIBRARY)(cgetMuc.o) \
		$(LIBRARY)(csampCS.o) \
		$(LIBRARY)(csampSoftMCS.o) \
		$(LIBRARY)(cMCSconstForMedia.o) \
		$(LIBRARY)(cfixMCSmodel.o)

