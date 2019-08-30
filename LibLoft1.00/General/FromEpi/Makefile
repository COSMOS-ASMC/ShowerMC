include  $(LIBLOFT)/site.config


EXTHDRS	      =

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = k16pGausLeg2.o \
		kisInPolyhed0.o \
		kbchop.o \
		kbetar.o \
		kbivarNorm.o \
		kbtest.o \
		kcsr1.o \
		kdgtest.o \
		kexpi.o \
		kklib.o \
		ksbeta.o \
		ksmpintbeta.o \
		ksort1.o \
		ksx2.o \
		psiim.o \
		kpolyg.o \
		kxplccl.o \
		kxplsq.o \
		kxplhorse.o \
		kseeinbox.o \
		kseeincyl.o \
		kxpLineAndTorus.o \
		kcquadEq.o \
		kquarticEq.o \
		kxplvsq.o \
		kinout2.o \
		kinout3.o \
		kgetNormalVec.o \
		kinout.o \
		kphiinside.o \
		krmLeadingB.o \
		ksplit2.o \
		kxplCylWall.o \
		kxplPipe.o \
		kgetAngleBtw2Planes.o \
		kProjAReaOfTriang.o \
		kxplPolyhed0.o \
		kxpLand3vpB.o \
		kCheckCoPlan.o \
		episoAngle.o \
		epcos2dir.o \
		afsep.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = k16pGausLeg2.f \
		kisInPolyhed0.f90 \
		kbchop.f \
		kbetar.f \
		kbivarNorm.f \
		kbtest.f \
		kcsr1.f \
		kdgtest.f \
		kexpi.f \
		kklib.f \
		ksbeta.f \
		ksmpintbeta.f \
		ksort1.f \
		ksx2.f \
		psiim.f \
		kpolyg.f \
		kxplccl.f \
		kxplsq.f \
		kxplhorse.f \
		kseeinbox.f \
		kseeincyl.f \
		kxpLineAndTorus.f90 \
		kcquadEq.f90 \
		kquarticEq.f90 \
		kxplvsq.f90  \
		kinout2.f90 \
		kinout3.f90 \
		kgetNormalVec.f90 \
		kinout.f90 \
		kphiinside.f90 \
		krmLeadingB.f90 \
		ksplit2.f90 \
		kxplCylWall.f90 \
		kxplPipe.f90 \
		kgetAngleBtw2Planes.f90 \
		kProjAReaOfTriang.f90  \
		kxplPolyhed0.f90 \
		kxpLand3vpB.f90 \
		kCheckCoPlan.f90 \
		episoAngle.f \
		epcos2dir.f90 \
		afsep.f

SYSHDRS	      =

all:
	make	Checkdir $(LIBRARY) SeeSubDir "TARGET=$@"

clean:
	@rm -f $(OBJS) core  *~ a.out 
	make SeeSubDir "TARGET=$@"

veryclean:
	@rm -f $(OBJS) core  *~ a.out temp*.f*
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
		cp -r  $(LIBLOFT)/Data  $(LIBLOFT)/Util \
		 $(LIBLOFT)/UserHook  $(LIBLOFT)/Scrpt  \
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
	                -a $$i != Data  -a $$i != Header \
	                -a $$i != Obso  -a $$i != TeX \
	                -a $$i != lib  -a $$i != UserHook \
	                -a $$i != Scrpt -a $$i != DistJob ] ; then \
		            $(MAKE) $(TARGET) \
			;fi) \
		; fi \
	done


$(LIBRARY):         \
		$(LIBRARY)(k16pGausLeg2.o) \
		$(LIBRARY)(kisInPolyhed0.o) \
	        $(LIBRARY)(kbchop.o) \
               $(LIBRARY)(kbetar.o) \
               $(LIBRARY)(kbivarNorm.o) \
               $(LIBRARY)(kbtest.o) \
               $(LIBRARY)(kcsr1.o) \
               $(LIBRARY)(kdgtest.o) \
               $(LIBRARY)(kexpi.o) \
               $(LIBRARY)(kklib.o) \
               $(LIBRARY)(ksbeta.o) \
               $(LIBRARY)(ksmpintbeta.o) \
               $(LIBRARY)(ksort1.o) \
               $(LIBRARY)(ksx2.o) \
               $(LIBRARY)(psiim.o)  \
               $(LIBRARY)(kpolyg.o)  \
               $(LIBRARY)(kxplccl.o)  \
               $(LIBRARY)(kxplsq.o)   \
		$(LIBRARY)(kxplhorse.o) \
		$(LIBRARY)(kseeinbox.o) \
		$(LIBRARY)(kseeincyl.o) \
		$(LIBRARY)(kxpLineAndTorus.o) \
		$(LIBRARY)(kcquadEq.o) \
		$(LIBRARY)(kquarticEq.o) \
		$(LIBRARY)(kxplvsq.o)  \
		$(LIBRARY)(kinout2.o) \
		$(LIBRARY)(kinout3.o) \
		$(LIBRARY)(kgetNormalVec.o) \
		$(LIBRARY)(kinout.o) \
		$(LIBRARY)(kphiinside.o) \
		$(LIBRARY)(krmLeadingB.o) \
		$(LIBRARY)(ksplit2.o)	\
		$(LIBRARY)(kxplCylWall.o) \
		$(LIBRARY)(kxplPipe.o) \
		$(LIBRARY)(kgetAngleBtw2Planes.o) \
		$(LIBRARY)(kProjAReaOfTriang.o) \
		$(LIBRARY)(kxplPolyhed0.o) \
		$(LIBRARY)(kxpLand3vpB.o) \
		$(LIBRARY)(kCheckCoPlan.o) \
		$(LIBRARY)(episoAngle.o) \
		$(LIBRARY)(epcos2dir.o) \
		$(LIBRARY)(afsep.o)







