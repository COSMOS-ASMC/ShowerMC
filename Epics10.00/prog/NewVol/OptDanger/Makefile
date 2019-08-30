include  $(EPICSTOP)/site.config


EXTHDRS	      = $(COSMOSINC)/Zglobalc.h \
		$(EPICSINC)/Zcnfig.h \
		$(EPICSINC)/Zep3Vec.h \
		$(EPICSINC)/ZepDirec.h \
		$(EPICSINC)/ZepMaxdef.h \
		$(EPICSINC)/ZepPos.h \
		$(EPICSINC)/ZepTrackp.h \
		$(EPICSINC)/Zepdebug.h



HDRS	      = Zciecone.h \
		Zelightg.h


MAKEFILE      = Makefile

OBJS	      = ep_octagon.o 

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = ep_octagon.f 
SYSHDRS	      =

all:
	make	Checkdir $(LIBRARY) SeeSubDir "TARGET=$@"

clean:
	@rm -f $(OBJS) core  *~ a.out 
	make SeeSubDir "TARGET=$@"

veryclean:
	@rm -f $(OBJS) core  *~ a.out temp*.f
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

print:
	@$(PRINT) $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"

tags:           $(HDRS) $(SRCS)
	@ctags $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"


Checkdir:
	@( cd $(EPICSTOP)/lib; if [ ! -d  $(ARCH) ]; then  mkdir $(ARCH); fi)


SeeSubDir:
	@for i in  * \
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



$(LIBRARY):         \
               $(LIBRARY)(ep_octagon.o) 
