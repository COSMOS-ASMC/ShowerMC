include  $(LIBLOFT)/site.config

EXTHDRS	      =

HDRS	      =

MAKEFILE      = Makefile


OBJS	      = epos-bas-lhc-v3400.o   \
		epos-con-lhc-v3400.o   \
		epos-dky-lhc-v3400.o  \
		epos-dro-lhc-v3400.o  \
		epos-ems-lhc-v3400.o  \
		epos-fra-lhc-v3400.o  \
		epos-hnb-lhc-v3400.o  \
		epos-ico-lhc-v3400.o  \
		epos-ids-lhc-v3400.o  \
		epos-int-lhc-v3400.o  \
		epos-jps-lhc-v3400.o  \
		eposm.o  \
		epos-omg-lhc-v3400.o  \
		epos-par-lhc-v3400.o  \
		epos-qsh-lhc-v3400.o  \
		epos-rsh-lhc-v3400.o  \
		epos-sem-lhc-v3400.o  \
		epos-tim-lhc-v3400.o  \
		eposu_no.o  \
		epos-uti-lhc-v3400.o  \
		epos-xan-lhc-v3400.o  \
		epos-xpr-lhc-v3400.o  \
		epos_j.o


PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = epos-bas-lhc-v3400.f  \
		epos-con-lhc-v3400.f  \
		epos-dky-lhc-v3400.f	\
		epos-dro-lhc-v3400.f  \
		epos-ems-lhc-v3400.f  \
		epos-fra-lhc-v3400.f  \
		epos-hnb-lhc-v3400.f  \
		epos-ico-lhc-v3400.f  \
		epos-ids-lhc-v3400.f  \
		epos-int-lhc-v3400.f  \
		epos-jps-lhc-v3400.f  \
		eposm.f  \
		epos-omg-lhc-v3400.f  \
		epos-par-lhc-v3400.f  \
		epos-qsh-lhc-v3400.f  \
		epos-rsh-lhc-v3400.f  \
		epos-sem-lhc-v3400.f  \
		epos-tim-lhc-v3400.f  \
		eposu_no.f  \
		epos-uti-lhc-v3400.f  \
		epos-xan-lhc-v3400.f  \
		epos-xpr-lhc-v3400.f  \
		epos_j.f 

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

$(LIBRARY):   \
		$(LIBRARY)(epos-bas-lhc-v3400.o)  \
		$(LIBRARY)(epos-con-lhc-v3400.o)  \
		$(LIBRARY)(epos-dky-lhc-v3400.o)  \
		$(LIBRARY)(epos-dro-lhc-v3400.o)  \
		$(LIBRARY)(epos-ems-lhc-v3400.o)  \
		$(LIBRARY)(epos-fra-lhc-v3400.o)  \
		$(LIBRARY)(epos-hnb-lhc-v3400.o)  \
		$(LIBRARY)(epos-ico-lhc-v3400.o)  \
		$(LIBRARY)(epos-ids-lhc-v3400.o)  \
		$(LIBRARY)(epos-int-lhc-v3400.o)  \
		$(LIBRARY)(epos-jps-lhc-v3400.o)  \
		$(LIBRARY)(eposm.o)  \
		$(LIBRARY)(epos-omg-lhc-v3400.o)  \
		$(LIBRARY)(epos-par-lhc-v3400.o)  \
		$(LIBRARY)(epos-qsh-lhc-v3400.o)  \
		$(LIBRARY)(epos-rsh-lhc-v3400.o)  \
		$(LIBRARY)(epos-sem-lhc-v3400.o)  \
		$(LIBRARY)(epos-tim-lhc-v3400.o)  \
		$(LIBRARY)(eposu_no.o)  \
		$(LIBRARY)(epos-uti-lhc-v3400.o)  \
		$(LIBRARY)(epos-xan-lhc-v3400.o)  \
		$(LIBRARY)(epos-xpr-lhc-v3400.o)  \
		$(LIBRARY)(epos_j.o)

