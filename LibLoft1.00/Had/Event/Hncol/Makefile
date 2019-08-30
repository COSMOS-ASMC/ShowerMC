include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcinippxc.h \
		$(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zevhnp.h \
		$(LIBLOFTINC)/Zevhnv.h \
		$(LIBLOFTINC)/Zmass.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/Zunionmap.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = cconsvChg.o \
		cgnlp.o \
		cgrap.o \
		chncol.o \
		cinipipx.o \
		cinippx.o \
		cs2lp.o \
		csampPolAng.o \
		csnchp.o \
		csptxy.o \
		cspwpt.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = cconsvChg.f \
		cgnlp.f \
		cgrap.f \
		chncol.f \
		cinipipx.f \
		cinippx.f \
		cs2lp.f \
		csampPolAng.f \
		csnchp.f \
		csptxy.f \
		cspwpt.f

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
$(LIBRARY)(cconsvChg.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cgnlp.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmass.h \
	$(LIBLOFTINC)/Zevhnv.h \
	$(LIBLOFTINC)/Zevhnp.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cgrap.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(chncol.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zevhnv.h \
	$(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cinipipx.o): $(LIBLOFTINC)/Zcinippxc.h
$(LIBRARY)(cinippx.o): $(LIBLOFTINC)/Zcinippxc.h
$(LIBRARY)(cs2lp.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmass.h \
	$(LIBLOFTINC)/Zevhnv.h \
	$(LIBLOFTINC)/Zcinippxc.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zevhnp.h
$(LIBRARY)(csnchp.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmass.h \
	$(LIBLOFTINC)/Zevhnv.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zevhnp.h
$(LIBRARY)(csptxy.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h

$(LIBRARY):         \
               $(LIBRARY)(cconsvChg.o) \
               $(LIBRARY)(cgnlp.o) \
               $(LIBRARY)(cgrap.o) \
               $(LIBRARY)(chncol.o) \
               $(LIBRARY)(cinipipx.o) \
               $(LIBRARY)(cinippx.o) \
               $(LIBRARY)(cs2lp.o) \
               $(LIBRARY)(csampPolAng.o) \
               $(LIBRARY)(csnchp.o) \
               $(LIBRARY)(csptxy.o) \
               $(LIBRARY)(cspwpt.o) 
