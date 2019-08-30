include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zcondc.h \
		$(LIBLOFTINC)/Zevhnp.h \
		$(LIBLOFTINC)/Zevhnv.h \
		$(LIBLOFTINC)/Zmass.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/Zunionmap.h \
	        $(LIBLOFTINC)/Zspecial.h \

HDRS	      = cumProb.h \
		cumProb14.h \
		cumProb16.h \
		cumProb2.h \
		cumProb40.h


MAKEFILE      = Makefile

OBJS	      = caveNoOfSucC.o \
		cfxTgtChg.o \
		chAcol.o \
		chAcolAdhoc.o \
		csampCollInA.o \
		csampFermiM.o \
		csmpColInA2.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = caveNoOfSucC.f \
		cfxTgtChg.f \
		chAcol.f \
		chAcolAdhoc.f \
		csampCollInA.f \
		csampFermiM.f \
		csmpColInA2.f

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
$(LIBRARY)(chAcol.o): $(LIBLOFTINC)/Zcondc.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zevhnp.h \
	$(LIBLOFTINC)/Zevhnv.h
$(LIBRARY)(chAcolAdhoc.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zevhnp.h
$(LIBRARY)(csampCollInA.o): $(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h cumProb.h
$(LIBRARY)(csampFermiM.o): $(LIBLOFTINC)/Zmass.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(csmpColInA2.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h cumProb2.h cumProb14.h \
	cumProb16.h cumProb40.h

$(LIBRARY):         \
               $(LIBRARY)(caveNoOfSucC.o) \
               $(LIBRARY)(cfxTgtChg.o) \
               $(LIBRARY)(chAcol.o) \
               $(LIBRARY)(chAcolAdhoc.o) \
               $(LIBRARY)(csampCollInA.o) \
               $(LIBRARY)(csampFermiM.o) \
               $(LIBRARY)(csmpColInA2.o) 
