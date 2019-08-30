include  $(COSMOSTOP)/site.config

LIBRARY		= $(LIBLOC)/libminuit.a

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = intrac.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = intrac.c

SYSHDRS	      =

all:
	make	Checkdir $(LIBRARY) SeeSubDir "TARGET=$@"
	@echo "*********************************************** *"
	@echo "* If Link process complains,                    *"
	@echo "* you may need to change the name of function   *"
	@echo "* in intrac.c in unix directory, say            *"
	@echo "* intrac_, _intrac, INTRAC, or ....          *"
	@echo "* depending on the system                    *"
	@echo "*                                            *"
	@echo "**********************************************"
	@echo "* I am sleeping  10 seconds                    *"
	@sleep 10
clean:
	@rm -f $(OBJS) core  *~ a.out   *.o
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
$(LIBRARY)(intrac..o):
$(LIBRARY):         \
               $(LIBRARY)(intrac.o)
