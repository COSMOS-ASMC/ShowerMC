include  $(LIBLOFT)/site.config


HDRS	      =


MAKEFILE      = Makefile

OBJS	      = kgetenv.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = kgetenv.c

SYSHDRS	      =

all:
	make	Checkdir $(LIBRARY) SeeSubDir "TARGET=$@"
#	@echo "**********************************************"
#	@echo "*                                            *"
#	@echo "* If current C compiler is 64 bit mode but   *"
#	@echo "* you want 32bit mode binary                 *"  
#	@echo "* you have to do the following by hand       *"
#	@echo "*   cd  Manager/Ccode                        *"
#	@echo "*   cc -c -m32 kgetenv.c                     *"
#	@echo "*   ar r $(LIBLOC)/libloft.a kgetenv.o     *"
#	@echo "*                                            *"
#	@echo "**********************************************"
#	@echo "  ***** I am sleeping 10 seconds *****        "
#	@sleep 10
 
clean:
	@rm -f $(OBJS) core  *~ a.out   *.o
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
$(LIBRARY)(kgetenv.o):
$(LIBRARY):         \
               $(LIBRARY)(kgetenv.o)
