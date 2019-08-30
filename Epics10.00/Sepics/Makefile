include  $(EPICSTOP)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zglobalc.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/Zunionmap.h \
		$(LIBLOFTINC)/Zep3Vec.h \
		$(LIBLOFTINC)/ZepDirec.h \
		$(LIBLOFTINC)/ZepManager.h \
		$(LIBLOFTINC)/ZepMaxdef.h \
		$(LIBLOFTINC)/ZepPos.h \
		$(LIBLOFTINC)/ZepTrackp.h \
		$(EPICSINC)/BlockData/epblksepi.h \
		$(EPICSINC)/Zcnfig.h \
		$(EPICSINC)/Zcnfig2.h \
		$(EPICSINC)/ZepTrack.h \
		$(EPICSINC)/ZsepManager.h \
		$(EPICSINC)/Zsparm.h \
		$(EPICSINC)/Zswk.h

HDRS	      =


MAKEFILE      = Makefile
#                 epshiftInciPos must be compiled first
OBJS	      = epShiftInciPos.o \
		epbdsep.o \
		sparmr.o \
		sepics.o 

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = epShiftInciPos.f90 \
		epbdsep.f \
		sparmr.f
#		sepics.f \

SYSHDRS	      =

all:
	make	Checkdir $(LIBRARY) SeeSubDir "TARGET=$@"

clean:
	@rm -f $(OBJS) core  *~ a.out 
	make SeeSubDir "TARGET=$@"

veryclean:
	@rm -f $(OBJS) core  *~ a.out temp*.f temp*.f90
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
		cp -r  $(EPICSTOP)/Data  $(EPICSTOP)/Util \
		 $(EPICSTOP)/UserHook  $(EPICSTOP)/Scrpt  \
		 $(EPICSTOP)/Doc        $(TEMPLATE)/
		cp $(EPICSTOP)/site.config $(TEMPLATE)/
		@echo "****" Installing header files into $(HEADER)
		cp -r $(EPICSTOP)/epics  $(HEADER)/
		@echo "****" Installation ended.
		@echo "****" Owner and group should be changed by hand

print:
	@$(PRINT) $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"

tags:           $(HDRS) $(SRCS)
	@ctags $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"


Checkdir:
	@( cd $(EPICSTOP)/lib; if [ ! -d  $(ARCH) ]; then  mkdir $(ARCH); fi)


SeeSubDir:
	@for i in Module * \
	;do   \
		if [ -d  $$i ] ; then \
			(cd $$i; \
			if [ -f Makefile  -a  $$i != Test \
			-a $$i != Util -a $$i != UC  \
	                -a $$i != Data  -a $$i != epics \
	                -a $$i != Obso  -a $$i != TeX \
	                -a $$i != lib  -a $$i != UserHook \
	                -a $$i != Scrpt -a $$i != DistJob ] ; then \
		            $(MAKE) $(TARGET) \
			;fi) \
		; fi \
	done


$(LIBRARY)(sepics.o): $(EPICSINC)/ZsepManager.h \
	$(LIBLOFTINC)/ZepDirec.h \
	$(LIBLOFTINC)/ZepTrackp.h \
	$(LIBLOFTINC)/ZepPos.h \
	$(LIBLOFTINC)/Zep3Vec.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/ZepManager.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zglobalc.h \
	$(EPICSINC)/Zcnfig.h \
	$(EPICSINC)/Zcnfig2.h \
	$(EPICSINC)/Zsparm.h \
	$(EPICSINC)/Zswk.h \
	$(EPICSINC)/ZepTrack.h 

$(LIBRARY)(sparmr.o):  \
	$(LIBLOFTINC)/ZepPos.h \
	$(LIBLOFTINC)/ZepDirec.h \
	$(LIBLOFTINC)/ZepManager.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(EPICSINC)/Zsparm.h \
	$(EPICSINC)/ZsepManager.h

$(LIBRARY):         \
               $(LIBRARY)(epShiftInciPos.o) \
               $(LIBRARY)(epbdsep.o) \
               $(LIBRARY)(sparmr.o) \
               $(LIBRARY)(sepics.o)
