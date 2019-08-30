include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/Zunionmap.h 

HDRS	      = 


MAKEFILE      = Makefile

OBJS	      = 	\
		ccos2gheCode.o \
		cghe2cos.o \
		cgheDecay.o \
		chAGheisha.o \
		cxsecGheisha.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = \
		ccos2gheCode.f \
		cghe2cos.f \
		cgheDecay.f \
		chAGheisha.f \
		cxsecGheisha.f 


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
$(LIBRARY)(ccos2gheCode.o): $(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(cghe2cos.o): $(LIBLOFTINC)/Zcode.h
$(LIBRARY)(cgheDecay.o): $(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h
$(LIBRARY)(chAGheisha.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h

$(LIBRARY)(cxsecGheisha.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h


$(LIBRARY):         \
		$(LIBRARY)(ccos2gheCode.o) \
	        $(LIBRARY)(cghe2cos.o) \
		$(LIBRARY)(cgheDecay.o) \
		$(LIBRARY)(chAGheisha.o) \
		$(LIBRARY)(cxsecGheisha.o)  

