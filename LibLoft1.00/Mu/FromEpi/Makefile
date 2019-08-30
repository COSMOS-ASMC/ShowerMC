include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zmass.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/Zunionmap.h \
		$(LIBLOFTINC)/ZbpPE.h \
		$(LIBLOFTINC)/ZbpSample.h \
		$(LIBLOFTINC)/ZbpTbl.h \
		$(LIBLOFTINC)/Zelement.h \
		$(LIBLOFTINC)/Zep3Vec.h \
		$(LIBLOFTINC)/ZepDirec.h \
		$(LIBLOFTINC)/ZepManager.h \
		$(LIBLOFTINC)/ZepMaxdef.h \
		$(LIBLOFTINC)/ZepPos.h \
		$(LIBLOFTINC)/ZepTrackp.h \
		$(LIBLOFTINC)/Zmedia.h \
		$(LIBLOFTINC)/Zmove.h \
		$(LIBLOFTINC)/Zmubpn.h \
		$(LIBLOFTINC)/Zstern.h \
		$(LIBLOFTINC)/Zurban.h \
		$(LIBLOFTINC)/Zmuint.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = epRdmuTab.o \
		epSetMu.o \
		epSetmuSTab.o \
		epmuBrsmp.o \
		epmuNsmp.o \
		epmuPrsmp.o \
		epmudEdx.o \
		epmuvmax2.o \
		eprdmuSmpTbl.o \
		epmuAuxP.o \
		epmudsdvdr.o \
		epmudsdv.o \
		epmuSetCnst.o \
		epmuvmax.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = epRdmuTab.f \
		epSetMu.f \
		epSetmuSTab.f \
		epmuBrsmp.f \
		epmuNsmp.f \
		epmuPrsmp.f \
		epmudEdx.f \
		epmuvmax2.f \
		eprdmuSmpTbl.f \
		epmuAuxP.f \
		epmudsdvdr.f \
		epmudsdv.f  \
		epmuSetCnst.f \
		epmuvmax.f  

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

install:	$(LIBRARY)
		@echo "****" Installing $(LIBRARY) into $(DEST)
		$(INSTALL) -f  $(DEST) $(LIBRARY)
		@echo "****" Installing Template files into $(TEMPLATE)
		cp -r  $(LIBLOFT)/Data  $(LIBLOFT)/Util \
		 $(LIBLOFT)/UserHook  $(LIBLOFT)/Scrpt  \
		 $(LIBLOFT)/Doc        $(TEMPLATE)/
		cp $(LIBLOFT)/site.config $(TEMPLATE)/
		@echo "****" Installing header files into $(HEADER)
		cp -r $(LIBLOFT)/epics  $(HEADER)/
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
	                -a $$i != Data  -a $$i != Header \
	                -a $$i != Obso  -a $$i != TeX \
	                -a $$i != lib  -a $$i != UserHook \
	                -a $$i != Scrpt -a $$i != DistJob ] ; then \
		            $(MAKE) $(TARGET) \
			;fi) \
		; fi \
	done

###
$(LIBRARY)(epRdmuTab.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/ZepTrackp.h
$(LIBRARY)(epSetMu.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epSetmuSTab.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epmuBrsmp.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epmuNsmp.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epmuPrsmp.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epmudEdx.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epmuvmax2.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(eprdmuSmpTbl.o): $(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/ZepManager.h \
	$(LIBLOFTINC)/ZepMaxdef.h

$(LIBRARY):         \
               $(LIBRARY)(epRdmuTab.o) \
               $(LIBRARY)(epSetMu.o) \
               $(LIBRARY)(epSetmuSTab.o) \
               $(LIBRARY)(epmuBrsmp.o) \
               $(LIBRARY)(epmuNsmp.o) \
               $(LIBRARY)(epmuPrsmp.o) \
               $(LIBRARY)(epmudEdx.o) \
               $(LIBRARY)(epmuvmax2.o) \
               $(LIBRARY)(eprdmuSmpTbl.o) \
               $(LIBRARY)(epmuAuxP.o) \
               $(LIBRARY)(epmudsdvdr.o) \
               $(LIBRARY)(epmudsdv.o) \
               $(LIBRARY)(epmuSetCnst.o) \
               $(LIBRARY)(epmuvmax.o) 
