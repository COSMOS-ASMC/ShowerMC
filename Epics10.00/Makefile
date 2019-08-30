include  $(EPICSTOP)/site.config


EXTHDRS	      =

HDRS	      =


MAKEFILE      = Makefile

OBJS	      =

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      =

SYSHDRS	      =

all:
	mkdir -p $(LIBLOC)	
	make	Checkdir $(LIBRARY) SeeSubDir "TARGET=$@"
	ranlib $(LIBRARY) 

clean:
	@rm -f $(OBJS) core *~
	make SeeSubDir "TARGET=$@"

veryclean:
	@rm -f $(OBJS) core *~ temp*.f temp*.f90
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
#		@echo "****" Installing Template files into $(TEMPLATE)
#		cp -r  $(EPICSTOP)/Data  $(EPICSTOP)/Util \
#		 $(EPICSTOP)/UserHook  $(EPICSTOP)/Scrpt \
#		 $(EPICSTOP)/Doc     $(TEMPLATE)/
#		cp $(EPICSTOP)/site.config $(TEMPLATE)/
#		@echo "****" Installing header files into $(HEADER)
#		cp -r $(EPICSTOP)/cosmos  $(HEADER)/
#		@echo "****" Installation ended.
#		@echo "****" Owner and group should be changed by hand

print:
	@$(PRINT) $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"

tags:           $(HDRS) $(SRCS)
	@ctags $(HDRS) $(SRCS)
	make SeeSubDir "TARGET=$@"


Checkdir:
	@( cd $(EPICSTOP)/lib; if [ ! -d  $(ARCH) ]; then  mkdir $(ARCH); fi)



SeeSubDir:
	@for i in Module  * \
	;do   \
		if [ -d  $$i ] ; then \
			(cd $$i; \
			if [ -f Makefile  -a  $$i != Test \
			-a $$i != Util -a $$i != UC  \
	                -a $$i != Ana  -a $$i != Ceren \
	                -a $$i != Data  -a $$i != epics \
	                -a $$i != Misc  -a $$i != MuUg \
	                -a $$i != Obso  -a $$i != TeX \
	                -a $$i != lib  -a $$i != UserHook \
	                -a $$i != Doc  -a $$i != Betterone \
	                -a $$i != Scrpt ] ; then \
		            $(MAKE) $(TARGET) \
			;fi) \
		; fi \
	done
	(cd Util/DrawNewVol;  $(MAKE) $(TARGET) )
$(LIBRARY):         
