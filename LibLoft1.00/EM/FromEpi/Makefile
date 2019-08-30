include  $(LIBLOFT)/site.config


EXTHDRS	      = $(LIBLOFTINC)/Zcode.h \
		$(LIBLOFTINC)/Zglobalc.h \
		$(LIBLOFTINC)/Zmass.h \
		$(LIBLOFTINC)/Zptcl.h \
		$(LIBLOFTINC)/Zunionmap.h \
		$(LIBLOFTINC)/BlockData/epblkemg.h \
		$(LIBLOFTINC)/ZbasicCnst.h \
		$(LIBLOFTINC)/ZbpPE.h \
		$(LIBLOFTINC)/ZbpSample.h \
		$(LIBLOFTINC)/ZbpTbl.h \
		$(LIBLOFTINC)/Zelement.h \
		$(LIBLOFTINC)/Zelmag.h \
		$(LIBLOFTINC)/Zep3Vec.h \
		$(LIBLOFTINC)/ZepDirec.h \
		$(LIBLOFTINC)/ZepManager.h \
		$(LIBLOFTINC)/ZepMaxdef.h \
		$(LIBLOFTINC)/ZepPos.h \
		$(LIBLOFTINC)/ZepTrack.h \
		$(LIBLOFTINC)/ZepTrackv.h \
		$(LIBLOFTINC)/Zepprm.h \
		$(LIBLOFTINC)/Zmedia.h \
		$(LIBLOFTINC)/Zmove.h \
		$(LIBLOFTINC)/Zmubpn.h \
		$(LIBLOFTINC)/Zstern.h \
		$(LIBLOFTINC)/Zurban.h

HDRS	      = ZdEdxSpec.h


MAKEFILE      = Makefile

OBJS	      = epSmpBremAng2BN.o \
		epSmpBremAngTsai.o \
		epSmpBremAngRelF.o \
		epBrCSamp.o \
		epBrHSamp.o \
		epBrLSamp.o \
		epBrSSamp.o \
		epBrSamp.o \
		epCompScrBr.o \
		epCompScrPr.o \
		epCoulombC.o \
		epExpot.o \
		epGetEffZA.o \
		epKnock.o \
		epPBA.o \
		epPrCSamp.o \
		epPrHSamp.o \
		epPrLSamp.o \
		epPrSamp.o \
		epPrTSamp.o \
		epReadMTbl.o \
		epReadTab.o \
		epReadXXsec.o \
		epSetPhotoE.o \
		epSetSTblCns.o \
		epStern.o \
		epUrban.o \
		epX0.o \
		epX0Old.o \
		epXrayp.o \
		epanihi.o \
		epbhabha.o \
		epcompton.o \
		epdedxNone.o \
		epdedxe.o \
		epdedxflhv.o \
		epdedxhvy.o \
		epmoller.o \
		epphotoE.o \
		eprdSmpTbl.o \
		epsync.o \
		epwtStern.o \
		epdEdxDenC.o \
		epMolfunc.o  \
		epMolBlogB.o  \
		epinteG.o  \
		epinteElec.o \
		epdoEMInte.o \
		epSampMolRedA.o \
		epSetTcut.o \
		epmulScat.o \
		epSampMol.o 

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = epSmpBremAng2BN.f90 \
		epSmpBremAngTsai.f90 \
		epSmpBremAngRelF.f90 \
		epBrCSamp.f \
		epBrHSamp.f \
		epBrLSamp.f \
		epBrSSamp.f \
		epBrSamp.f \
		epCompScrBr.f \
		epCompScrPr.f \
		epCoulombC.f \
		epExpot.f \
		epGetEffZA.f \
		epKnock.f \
		epPBA.f \
		epPrCSamp.f \
		epPrHSamp.f \
		epPrLSamp.f \
		epPrSamp.f \
		epPrTSamp.f \
		epReadMTbl.f \
		epReadTab.f \
		epReadXXsec.f \
		epSetPhotoE.f \
		epSetSTblCns.f \
		epStern.f \
		epUrban.f \
		epX0.f \
		epX0Old.f \
		epXrayp.f \
		epanihi.f \
		epbhabha.f \
		epcompton.f \
		epdedxNone.f \
		epdedxe.f \
		epdedxflhv.f \
		epdedxhvy.f \
		epmoller.f \
		epphotoE.f \
		eprdSmpTbl.f \
		epsync.f \
		epwtStern.f \
		epdEdxDenC.f \
		epMolfunc.f90 \
		epMolBlogB.f90 \
		epinteG.f \
		epinteElec.f \
		epdoEMInte.f \
		epSampMolRedA.f90 \
		epSetTcut.f  \
		epmulScat.f  \
		epSampMol.f
SYSHDRS	      =

all:
	make	Checkdir $(LIBRARY) SeeSubDir "TARGET=$@"

clean:
	@rm -f $(OBJS) core  *~ a.out 
	make SeeSubDir "TARGET=$@"

veryclean:
	@rm -f $(OBJS) core  *~ a.out temp*.f temp*.f90
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

###
$(LIBRARY)(epSmpBremAng2BN.o): $(LIBLOFTINC)/Zmass.h

$(LIBRARY)(epSmpBremAngTsai.o): $(LIBLOFTINC)/Zmass.h




$(LIBRARY)(epBrCSamp.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/ZbasicCnst.h \
	$(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epBrHSamp.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epBrLSamp.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epBrSSamp.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epBrSamp.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epCompScrBr.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/ZbasicCnst.h
$(LIBRARY)(epCompScrPr.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/ZbasicCnst.h
$(LIBRARY)(epExpot.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epGetEffZA.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/ZbasicCnst.h \
	$(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epKnock.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epPrCSamp.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/ZbasicCnst.h \
	$(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epPrHSamp.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zglobalc.h
$(LIBRARY)(epPrLSamp.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epPrSamp.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epPrTSamp.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h \
	$(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/ZbasicCnst.h
$(LIBRARY)(epReadMTbl.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epReadTab.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epReadXXsec.o): $(LIBLOFTINC)/ZepManager.h \
	$(LIBLOFTINC)/ZepMaxdef.h
$(LIBRARY)(epSetPhotoE.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epSetSTblCns.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epStern.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epUrban.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zurban.h ZdEdxSpec.h \
	$(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epX0.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epXrayp.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/ZbasicCnst.h \
	$(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h
$(LIBRARY)(epanihi.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epbhabha.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epcompton.o): $(LIBLOFTINC)/Zglobalc.h \
	$(LIBLOFTINC)/ZbasicCnst.h \
	$(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epdedxNone.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmass.h ZdEdxSpec.h
$(LIBRARY)(epdedxe.o): $(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zmass.h ZdEdxSpec.h
$(LIBRARY)(epdedxflhv.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epdedxhvy.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epmoller.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(epphotoE.o): $(LIBLOFTINC)/Zmedia.h \
	$(LIBLOFTINC)/ZepMaxdef.h \
	$(LIBLOFTINC)/Zelement.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/Zstern.h \
	$(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpPE.h \
	$(LIBLOFTINC)/Zurban.h \
	$(LIBLOFTINC)/Zmubpn.h \
	$(LIBLOFTINC)/Zmass.h
$(LIBRARY)(eprdSmpTbl.o): $(LIBLOFTINC)/ZbpSample.h \
	$(LIBLOFTINC)/ZbpTbl.h \
	$(LIBLOFTINC)/ZepManager.h \
	$(LIBLOFTINC)/ZepMaxdef.h
$(LIBRARY)(epsync.o): $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zunionmap.h \
	$(LIBLOFTINC)/Zep3Vec.h
$(LIBRARY)(epwtStern.o): $(LIBLOFTINC)/Zstern.h

$(LIBRARY)(epdEdxDenC.o): $(LIBLOFTINC)/Zstern.h \
		$(LIBLOFTINC)/Zmedia.h
$(LIBRARY)(epSampMol.o): $(LIBLOFTINC)/Zmedia.h  \
                $(LIBLOFTINC)/Zcode.h  \
                $(LIBLOFTINC)/Zmass.h \
		$(LIBLOFTINC)/Zglobalc.h
$(LIBRARY)(epmulScat.o): $(LIBLOFTINC)/Zmedia.h  \
                $(LIBLOFTINC)/Zcode.h  \
                $(LIBLOFTINC)/Zmass.h \
		$(LIBLOFTINC)/Zglobalc.h


$(LIBRARY):         \
               $(LIBRARY)(epSmpBremAng2BN.o) \
               $(LIBRARY)(epSmpBremAngTsai.o) \
               $(LIBRARY)(epSmpBremAngRelF.o) \
               $(LIBRARY)(epBrCSamp.o) \
               $(LIBRARY)(epBrHSamp.o) \
               $(LIBRARY)(epBrLSamp.o) \
               $(LIBRARY)(epBrSSamp.o) \
               $(LIBRARY)(epBrSamp.o) \
               $(LIBRARY)(epCompScrBr.o) \
               $(LIBRARY)(epCompScrPr.o) \
               $(LIBRARY)(epCoulombC.o) \
               $(LIBRARY)(epExpot.o) \
               $(LIBRARY)(epGetEffZA.o) \
               $(LIBRARY)(epKnock.o) \
               $(LIBRARY)(epPBA.o) \
               $(LIBRARY)(epPrCSamp.o) \
               $(LIBRARY)(epPrHSamp.o) \
               $(LIBRARY)(epPrLSamp.o) \
               $(LIBRARY)(epPrSamp.o) \
               $(LIBRARY)(epPrTSamp.o) \
               $(LIBRARY)(epReadMTbl.o) \
               $(LIBRARY)(epReadTab.o) \
               $(LIBRARY)(epReadXXsec.o) \
               $(LIBRARY)(epSetPhotoE.o) \
               $(LIBRARY)(epSetSTblCns.o) \
               $(LIBRARY)(epStern.o) \
               $(LIBRARY)(epUrban.o) \
               $(LIBRARY)(epX0.o) \
               $(LIBRARY)(epX0Old.o) \
               $(LIBRARY)(epXrayp.o) \
               $(LIBRARY)(epanihi.o) \
               $(LIBRARY)(epbhabha.o) \
               $(LIBRARY)(epcompton.o) \
               $(LIBRARY)(epdedxNone.o) \
               $(LIBRARY)(epdedxe.o) \
               $(LIBRARY)(epdedxflhv.o) \
               $(LIBRARY)(epdedxhvy.o) \
               $(LIBRARY)(epmoller.o) \
               $(LIBRARY)(epphotoE.o) \
               $(LIBRARY)(eprdSmpTbl.o) \
               $(LIBRARY)(epsync.o) \
               $(LIBRARY)(epwtStern.o) \
               $(LIBRARY)(epdEdxDenC.o) \
               $(LIBRARY)(epMolfunc.o) \
               $(LIBRARY)(epMolBlogB.o) \
               $(LIBRARY)(epinteG.o) \
               $(LIBRARY)(epinteElec.o) \
               $(LIBRARY)(epdoEMInte.o)  \
               $(LIBRARY)(epSampMolRedA.o) \
               $(LIBRARY)(epSetTcut.o) \
		$(LIBRARY)(epmulScat.o) \
		$(LIBRARY)(epSampMol.o)



