include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zkcele.h

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = k16pGaussLeg.o \
		k4ptdi.o \
		kbetar.o \
		kbinChop.o \
		kbinChopWP.o \
		kSbinChopWP.o \
		kbinom.o \
		kceles.o \
		kcombi.o \
		kcossn.o \
		kcsplCoef.o \
		kNewtonRaphson.o \
		kSNewtonRaphson.o \
		kScsplCoef.o \
		kcsplDif.o \
		kScsplDif.o \
		kcsplInteg.o \
		kScsplInteg.o \
		kcsplIntp.o \
		kcsr1idx.o \
		kScsplIntp.o \
		kdexpIntF.o \
		kdexpIntFWP.o \
		kdexpIntFb.o \
		kdir2deg.o \
		kdmachmnmx.o \
		kdwhereis.o \
		kerf.o \
		kerfc.o \
		kexpot.o \
		kfrge.o \
		kgamma.o \
		kgauss.o \
		kgsub.o \
		kintp3.o \
		klena.o \
		kmanbit.o \
		knbino.o \
		kpmnx.o \
		kpoisn.o \
		kpolintp.o \
		kpolintp2.o \
		kpolintpS.o \
		kpolintp2S.o \
		kqsort.o \
		kqsortc.o \
		kqsortd.o \
		kqsorth.o \
		kqsorti.o \
		kqsortr.o \
		kroundrv.o \
		ksampLin.o \
		ksampPEang.o \
		ksampPw.o \
		ksampRSA.o \
		ksbwig.o \
		kseblk.o \
		ksgamd.o \
		ksmpConExp.o \
		ksortinv.o \
		ksplandau.o \
		ksplit.o \
		ksupblank.o \
		ktrpzIntT2.o \
		ktrpzIntegT.o \
		kzeta2a.o \
		ksmooth.o \
		ksampAF.o \
		kcountFields.o \
		kSRunge-Kutta4.o \
		kvec_prod3.o \
		ksampForBPA.o \
		kcubicEq.o \
		kpolygC.o \
		kexpiC.o \
		kquadEq.o \
		rnd.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = k16pGaussLeg.f \
		k4ptdi.f \
		kbetar.f \
		kbinChop.f \
		kbinChopWP.f90 \
		kSbinChopWP.f90 \
		kbinom.f \
		kceles.f \
		kcombi.f \
		kcossn.f \
		kcsplCoef.f \
		kNewtonRaphson.f90 \
		kSNewtonRaphson.f90 \
		kScsplCoef.f \
		kcsplDif.f \
		kScsplDif.f \
		kcsplInteg.f \
		kScsplInteg.f \
		kcsplIntp.f \
		kScsplIntp.f \
		kdexpIntF.f \
		kdexpIntFWP.f90 \
		kdexpIntFb.f \
		kdir2deg.f \
		kdmachmnmx.f \
		kdwhereis.f \
		kerf.f \
		kerfc.f \
		kexpot.f \
		kfrge.f \
		kgamma.f \
		kgauss.f \
		kgsub.f \
		kintp3.f \
		klena.f \
		kmanbit.f \
		knbino.f \
		kpmnx.f \
		kpoisn.f \
		kpolintp.f \
		kpolintp2.f \
		kpolintpS.f \
		kpolintp2S.f \
		kqsort.f \
		kqsortc.f \
		kqsortd.f \
		kqsorth.f \
		kqsorti.f \
		kqsortr.f \
		kroundrv.f \
		ksampLin.f \
		ksampPEang.f \
		ksampPw.f \
		ksampRSA.f \
		ksbwig.f \
		kseblk.f \
		ksgamd.f \
		ksmpConExp.f \
		ksortinv.f \
		ksplandau.f \
		ksplit.f \
		ksupblank.f \
		ktrpzIntT2.f \
		ktrpzIntegT.f \
		kzeta2a.f \
		ksmooth.f \
		ksampAF.f \
		kcountFields.f90 \
		kSRunge-Kutta4.f90 \
		kvec_prod3.f90 \
		ksampForBPA.f90 \
		kcubicEq.f90 \
		kpolygC.f \
		kexpiC.f \
		kcsr1idx.f \
		kquadEq.f90 \
		rnd.f

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
$(LIBRARY)(kceles.o): $(LIBLOFTINC)/Zkcele.h

$(LIBRARY):         \
               $(LIBRARY)(k16pGaussLeg.o) \
               $(LIBRARY)(k4ptdi.o) \
               $(LIBRARY)(kbetar.o) \
               $(LIBRARY)(kbinChop.o) \
               $(LIBRARY)(kbinChopWP.o) \
               $(LIBRARY)(kSbinChopWP.o) \
               $(LIBRARY)(kbinom.o) \
               $(LIBRARY)(kceles.o) \
               $(LIBRARY)(kcombi.o) \
               $(LIBRARY)(kcossn.o) \
               $(LIBRARY)(kcsplCoef.o) \
               $(LIBRARY)(kNewtonRaphson.o) \
               $(LIBRARY)(kSNewtonRaphson.o) \
               $(LIBRARY)(kScsplCoef.o) \
               $(LIBRARY)(kcsplDif.o) \
               $(LIBRARY)(kScsplDif.o) \
               $(LIBRARY)(kcsplInteg.o) \
               $(LIBRARY)(kScsplInteg.o) \
               $(LIBRARY)(kcsplIntp.o) \
               $(LIBRARY)(kcsr1idx.o) \
               $(LIBRARY)(kScsplIntp.o) \
               $(LIBRARY)(kdexpIntF.o) \
               $(LIBRARY)(kdexpIntFWP.o) \
               $(LIBRARY)(kdexpIntFb.o) \
               $(LIBRARY)(kdir2deg.o) \
               $(LIBRARY)(kdmachmnmx.o) \
               $(LIBRARY)(kdwhereis.o) \
               $(LIBRARY)(kerf.o) \
               $(LIBRARY)(kerfc.o) \
               $(LIBRARY)(kexpot.o) \
               $(LIBRARY)(kfrge.o) \
               $(LIBRARY)(kgamma.o) \
               $(LIBRARY)(kgauss.o) \
               $(LIBRARY)(kgsub.o) \
               $(LIBRARY)(kintp3.o) \
               $(LIBRARY)(klena.o) \
               $(LIBRARY)(kmanbit.o) \
               $(LIBRARY)(knbino.o) \
               $(LIBRARY)(kpmnx.o) \
               $(LIBRARY)(kpoisn.o) \
               $(LIBRARY)(kpolintp.o) \
               $(LIBRARY)(kpolintp2.o) \
               $(LIBRARY)(kpolintpS.o) \
               $(LIBRARY)(kpolintp2S.o) \
               $(LIBRARY)(kqsort.o) \
               $(LIBRARY)(kqsortc.o) \
               $(LIBRARY)(kqsortd.o) \
               $(LIBRARY)(kqsorth.o) \
               $(LIBRARY)(kqsorti.o) \
               $(LIBRARY)(kqsortr.o) \
               $(LIBRARY)(kroundrv.o) \
               $(LIBRARY)(ksampLin.o) \
               $(LIBRARY)(ksampPEang.o) \
               $(LIBRARY)(ksampPw.o) \
               $(LIBRARY)(ksampRSA.o) \
               $(LIBRARY)(ksbwig.o) \
               $(LIBRARY)(kseblk.o) \
               $(LIBRARY)(ksgamd.o) \
               $(LIBRARY)(ksmpConExp.o) \
               $(LIBRARY)(ksortinv.o) \
               $(LIBRARY)(ksplandau.o) \
               $(LIBRARY)(ksplit.o) \
               $(LIBRARY)(ksupblank.o) \
               $(LIBRARY)(ktrpzIntT2.o) \
               $(LIBRARY)(ktrpzIntegT.o) \
               $(LIBRARY)(kzeta2a.o) \
               $(LIBRARY)(ksmooth.o) \
               $(LIBRARY)(ksampAF.o) \
               $(LIBRARY)(kcountFields.o) \
               $(LIBRARY)(kSRunge-Kutta4.o) \
               $(LIBRARY)(kvec_prod3.o) \
               $(LIBRARY)(ksampForBPA.o) \
               $(LIBRARY)(kcubicEq.o) \
               $(LIBRARY)(kpolygC.o) \
               $(LIBRARY)(kexpiC.o) \
               $(LIBRARY)(kquadEq.o) \
		$(LIBRARY)(rnd.o)
