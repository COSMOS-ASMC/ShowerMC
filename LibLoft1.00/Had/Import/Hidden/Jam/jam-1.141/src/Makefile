include  $(LIBLOFT)/site.config


EXTHDRS	      =

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = gfortran.o \
		jam.o \
		jamana.o \
		jamcross.o \
		jambuu.o \
		jamcoll.o \
		jamdec.o \
		jamdummy.o \
		jamhard.o \
		jamhij.o \
		jampdf.o \
		jamsoft.o \
		jampyjet.o \
		jampysigh.o \
		jampythia.o  \
		jamdat.o


PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = gfortran.f \
		jam.f \
		jamana.f \
		jambuu.f \
		jamcoll.f \
		jamcross.f \
		jamdec.f \
		jamdummy.f \
		jamhard.f \
		jamhij.f \
		jampdf.f \
		jamsoft.f \
		jampyjet.f \
		jampysigh.f \
		jampythia.f  \
		jamdat.f 


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

$(LIBRARY):         \
		$(LIBRARY)(gfortran.o) \
		$(LIBRARY)(jam.o) \
		$(LIBRARY)(jamana.o) \
		$(LIBRARY)(jambuu.o) \
		$(LIBRARY)(jamcoll.o) \
		$(LIBRARY)(jamcross.o) \
		$(LIBRARY)(jamdec.o) \
		$(LIBRARY)(jamdummy.o) \
		$(LIBRARY)(jamhard.o) \
		$(LIBRARY)(jamhij.o) \
		$(LIBRARY)(jampdf.o) \
		$(LIBRARY)(jamsoft.o) \
		$(LIBRARY)(jampyjet.o) \
		$(LIBRARY)(jampysigh.o) \
		$(LIBRARY)(jampythia.o)  \
		$(LIBRARY)(jamdat.o)
