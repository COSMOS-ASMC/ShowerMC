include  $(COSMOSTOP)/site.config
LIBRARY	      = $(LIBLOC)/libminuit.a

HDRS	      =

MAKEFILE      = Makefile

OBJS	      =  minuit.o \
		 mnamin.o \
		 mnbins.o \
		 mncalf.o \
		 mncler.o \
		 mncntr.o \
		 mncomd.o \
		 mncont.o \
		 mncrck.o \
		 mncros.o \
		 mncuve.o \
		 mnderi.o \
		 mndxdi.o \
		 mneig.o \
		 mnemat.o \
		 mnerrs.o \
		 mneval.o \
		 mnexcm.o \
		 mnexin.o \
		 mnfixp.o \
		 mnfree.o \
		 mngrad.o \
		 mnhelp.o \
		 mnhes1.o \
		 mnhess.o \
		 mnimpr.o \
		 mninex.o \
		 mninit.o \
		 mninpu.o \
		 mnintr.o \
		 mnlims.o \
		 mnline.o \
		 mnmatu.o \
		 mnmigr.o \
		 mnmnos.o \
		 mnmnot.o \
		 mnparm.o \
		 mnpars.o \
		 mnpfit.o \
		 mnpint.o \
		 mnplot.o \
		 mnpout.o \
		 mnprin.o \
		 mnpsdf.o \
		 mnrazz.o \
		 mnread.o \
		 mnrn15.o \
		 mnrset.o \
		 mnsave.o \
		 mnscan.o \
		 mnseek.o \
		 mnset.o \
		 mnseti.o \
		 mnsimp.o \
		 mnstat.o \
		 mnstin.o \
		 mntiny.o \
		 mnunpt.o \
		 mnvers.o \
		 mnvert.o \
		 mnwarn.o \
		 mnwerr.o \
		 stand.o 


PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = minuit.f \
		mnamin.f \
		mnbins.f \
		mncalf.f \
		mncler.f \
		mncntr.f \
		mncomd.f \
		mncont.f \
		mncrck.f \
		mncros.f \
		mncuve.f \
		mnderi.f \
		mndxdi.f \
		mneig.f \
		mnemat.f \
		mnerrs.f \
		mneval.f \
		mnexcm.f \
		mnexin.f \
		mnfixp.f \
		mnfree.f \
		mngrad.f \
		mnhelp.f \
		mnhes1.f \
		mnhess.f \
		mnimpr.f \
		mninex.f \
		mninit.f \
		mninpu.f \
		mnintr.f \
		mnlims.f \
		mnline.f \
		mnmatu.f \
		mnmigr.f \
		mnmnos.f \
		mnmnot.f \
		mnparm.f \
		mnpars.f \
		mnpfit.f \
		mnpint.f \
		mnplot.f \
		mnpout.f \
		mnprin.f \
		mnpsdf.f \
		mnrazz.f \
		mnread.f \
		mnrn15.f \
		mnrset.f \
		mnsave.f \
		mnscan.f \
		mnseek.f \
		mnset.f \
		mnseti.f \
		mnsimp.f \
		mnstat.f \
		mnstin.f \
		mntiny.f \
		mnunpt.f \
		mnvers.f \
		mnvert.f \
		mnwarn.f \
		mnwerr.f \
		stand.f 

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

$(LIBRARY):         \
		$(LIBRARY)(minuit.o)  \
		$(LIBRARY)(mnamin.o)  \
		$(LIBRARY)(mnbins.o)  \
		$(LIBRARY)(mncalf.o)  \
		$(LIBRARY)(mncler.o)  \
		$(LIBRARY)(mncntr.o)  \
		$(LIBRARY)(mncomd.o)  \
		$(LIBRARY)(mncont.o)  \
		$(LIBRARY)(mncrck.o)  \
		$(LIBRARY)(mncros.o)  \
		$(LIBRARY)(mncuve.o)  \
		$(LIBRARY)(mnderi.o)  \
		$(LIBRARY)(mndxdi.o)  \
		$(LIBRARY)(mneig.o)  \
		$(LIBRARY)(mnemat.o)  \
		$(LIBRARY)(mnerrs.o)  \
		$(LIBRARY)(mneval.o)  \
		$(LIBRARY)(mnexcm.o)  \
		$(LIBRARY)(mnexin.o)  \
		$(LIBRARY)(mnfixp.o)  \
		$(LIBRARY)(mnfree.o)  \
		$(LIBRARY)(mngrad.o)  \
		$(LIBRARY)(mnhelp.o)  \
		$(LIBRARY)(mnhes1.o)  \
		$(LIBRARY)(mnhess.o)  \
		$(LIBRARY)(mnimpr.o)  \
		$(LIBRARY)(mninex.o)  \
		$(LIBRARY)(mninit.o)  \
		$(LIBRARY)(mninpu.o)  \
		$(LIBRARY)(mnintr.o)  \
		$(LIBRARY)(mnlims.o)  \
		$(LIBRARY)(mnline.o)  \
		$(LIBRARY)(mnmatu.o)  \
		$(LIBRARY)(mnmigr.o)  \
		$(LIBRARY)(mnmnos.o)  \
		$(LIBRARY)(mnmnot.o)  \
		$(LIBRARY)(mnparm.o)  \
		$(LIBRARY)(mnpars.o)  \
		$(LIBRARY)(mnpfit.o)  \
		$(LIBRARY)(mnpint.o)  \
		$(LIBRARY)(mnplot.o)  \
		$(LIBRARY)(mnpout.o)  \
		$(LIBRARY)(mnprin.o)  \
		$(LIBRARY)(mnpsdf.o)  \
		$(LIBRARY)(mnrazz.o)  \
		$(LIBRARY)(mnread.o)  \
		$(LIBRARY)(mnrn15.o)  \
		$(LIBRARY)(mnrset.o)  \
		$(LIBRARY)(mnsave.o)  \
		$(LIBRARY)(mnscan.o)  \
		$(LIBRARY)(mnseek.o)  \
		$(LIBRARY)(mnset.o)  \
		$(LIBRARY)(mnseti.o)  \
		$(LIBRARY)(mnsimp.o)  \
		$(LIBRARY)(mnstat.o)  \
		$(LIBRARY)(mnstin.o)  \
		$(LIBRARY)(mntiny.o)  \
		$(LIBRARY)(mnunpt.o)  \
		$(LIBRARY)(mnvers.o)  \
		$(LIBRARY)(mnvert.o)  \
		$(LIBRARY)(mnwarn.o)  \
		$(LIBRARY)(mnwerr.o)  \
		$(LIBRARY)(stand.o)  
