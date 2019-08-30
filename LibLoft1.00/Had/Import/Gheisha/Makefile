include  $(LIBLOFT)/site.config


EXTHDRS	      =

HDRS	      =


MAKEFILE      = Makefile

OBJS	      = acfsin.o \
		add.o \
		ang.o \
		atomas.o \
		captur.o \
		casal0.o \
		casaom.o \
		casasm.o \
		casasp.o \
		casax0.o \
		casaxm.o \
		casfrg.o \
		cask0.o \
		cask0b.o \
		caskm.o \
		caskp.o \
		casl0.o \
		casn.o \
		casnb.o \
		casom.o \
		casp.o \
		caspb.o \
		caspim.o \
		caspip.o \
		cassm.o \
		cassp.o \
		casx0.o \
		casxm.o \
		chanwx.o \
		cinema.o \
		compo.o \
		coranh.o \
		coscat.o \
		defs.o \
		defs1.o \
		dlpnuc.o \
		dlpsor.o \
		dotes2.o \
		dotnuc.o \
		dpdnuc.o \
		exnu.o \
		fctcos.o \
		fermi.o \
		fissio.o \
		gdpdk.o \
		genxpt.o \
		gheini.o \
		gheish.o \
		ghstop.o \
		ghtype.o \
		gmate.o \
		gnslwd.o \
		gpart.o \
		gpdk.o \
		gpghei.o \
		gpions.o \
		gspart.o \
		gutil.o \
		hadini.o \
		hadnuc.o \
		kmabs.o \
		lengtx.o \
		lor.o \
		minv.o \
		nbanh.o \
		nfac.o \
		normal.o \
		nucrec.o \
		pbanh.o \
		pcsdat.o \
		pdgtyp.o \
		phasp.o \
		phpnuc.o \
		pimabs.o \
		pmltpc.o \
		poisso.o \
		ranres.o \
		rotes2.o \
		rtmi.o \
		selfab.o \
		setcur.o \
		settrk.o \
		sginit.o \
		steep.o \
		steeq.o \
		steer.o \
		stees.o \
		stpair.o \
		sub.o \
		tdelay.o \
		trac.o \
		twob.o \
		twoclu.o

PRINT	      = pr

SHELL	      = /bin/sh

SRCS	      = acfsin.f \
		add.f \
		ang.f \
		atomas.f \
		captur.f \
		casal0.f \
		casaom.f \
		casasm.f \
		casasp.f \
		casax0.f \
		casaxm.f \
		casfrg.f \
		cask0.f \
		cask0b.f \
		caskm.f \
		caskp.f \
		casl0.f \
		casn.f \
		casnb.f \
		casom.f \
		casp.f \
		caspb.f \
		caspim.f \
		caspip.f \
		cassm.f \
		cassp.f \
		casx0.f \
		casxm.f \
		chanwx.f \
		cinema.f \
		compo.f \
		coranh.f \
		coscat.f \
		defs.f \
		defs1.f \
		dlpnuc.f \
		dlpsor.f \
		dotes2.f \
		dotnuc.f \
		dpdnuc.f \
		exnu.f \
		fctcos.f \
		fermi.f \
		fissio.f \
		gdpdk.f \
		genxpt.f \
		gheini.f \
		gheish.f \
		ghstop.f \
		ghtype.f \
		gmate.f \
		gnslwd.f \
		gpart.f \
		gpdk.f \
		gpghei.f \
		gpions.f \
		gspart.f \
		gutil.f \
		hadini.f \
		hadnuc.f \
		kmabs.f \
		lengtx.f \
		lor.f \
		minv.f \
		nbanh.f \
		nfac.f \
		normal.f \
		nucrec.f \
		pbanh.f \
		pcsdat.f \
		pdgtyp.f \
		phasp.f \
		phpnuc.f \
		pimabs.f \
		pmltpc.f \
		poisso.f \
		ranres.f \
		rotes2.f \
		rtmi.f \
		selfab.f \
		setcur.f \
		settrk.f \
		sginit.f \
		steep.f \
		steeq.f \
		steer.f \
		stees.f \
		stpair.f \
		sub.f \
		tdelay.f \
		trac.f \
		twob.f \
		twoclu.f

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
               $(LIBRARY)(acfsin.o) \
               $(LIBRARY)(add.o) \
               $(LIBRARY)(ang.o) \
               $(LIBRARY)(atomas.o) \
               $(LIBRARY)(captur.o) \
               $(LIBRARY)(casal0.o) \
               $(LIBRARY)(casaom.o) \
               $(LIBRARY)(casasm.o) \
               $(LIBRARY)(casasp.o) \
               $(LIBRARY)(casax0.o) \
               $(LIBRARY)(casaxm.o) \
               $(LIBRARY)(casfrg.o) \
               $(LIBRARY)(cask0.o) \
               $(LIBRARY)(cask0b.o) \
               $(LIBRARY)(caskm.o) \
               $(LIBRARY)(caskp.o) \
               $(LIBRARY)(casl0.o) \
               $(LIBRARY)(casn.o) \
               $(LIBRARY)(casnb.o) \
               $(LIBRARY)(casom.o) \
               $(LIBRARY)(casp.o) \
               $(LIBRARY)(caspb.o) \
               $(LIBRARY)(caspim.o) \
               $(LIBRARY)(caspip.o) \
               $(LIBRARY)(cassm.o) \
               $(LIBRARY)(cassp.o) \
               $(LIBRARY)(casx0.o) \
               $(LIBRARY)(casxm.o) \
               $(LIBRARY)(chanwx.o) \
               $(LIBRARY)(cinema.o) \
               $(LIBRARY)(compo.o) \
               $(LIBRARY)(coranh.o) \
               $(LIBRARY)(coscat.o) \
               $(LIBRARY)(defs.o) \
               $(LIBRARY)(defs1.o) \
               $(LIBRARY)(dlpnuc.o) \
               $(LIBRARY)(dlpsor.o) \
               $(LIBRARY)(dotes2.o) \
               $(LIBRARY)(dotnuc.o) \
               $(LIBRARY)(dpdnuc.o) \
               $(LIBRARY)(exnu.o) \
               $(LIBRARY)(fctcos.o) \
               $(LIBRARY)(fermi.o) \
               $(LIBRARY)(fissio.o) \
               $(LIBRARY)(gdpdk.o) \
               $(LIBRARY)(genxpt.o) \
               $(LIBRARY)(gheini.o) \
               $(LIBRARY)(gheish.o) \
               $(LIBRARY)(ghstop.o) \
               $(LIBRARY)(ghtype.o) \
               $(LIBRARY)(gmate.o) \
               $(LIBRARY)(gnslwd.o) \
               $(LIBRARY)(gpart.o) \
               $(LIBRARY)(gpdk.o) \
               $(LIBRARY)(gpghei.o) \
               $(LIBRARY)(gpions.o) \
               $(LIBRARY)(gspart.o) \
               $(LIBRARY)(gutil.o) \
               $(LIBRARY)(hadini.o) \
               $(LIBRARY)(hadnuc.o) \
               $(LIBRARY)(kmabs.o) \
               $(LIBRARY)(lengtx.o) \
               $(LIBRARY)(lor.o) \
               $(LIBRARY)(minv.o) \
               $(LIBRARY)(nbanh.o) \
               $(LIBRARY)(nfac.o) \
               $(LIBRARY)(normal.o) \
               $(LIBRARY)(nucrec.o) \
               $(LIBRARY)(pbanh.o) \
               $(LIBRARY)(pcsdat.o) \
               $(LIBRARY)(pdgtyp.o) \
               $(LIBRARY)(phasp.o) \
               $(LIBRARY)(phpnuc.o) \
               $(LIBRARY)(pimabs.o) \
               $(LIBRARY)(pmltpc.o) \
               $(LIBRARY)(poisso.o) \
               $(LIBRARY)(ranres.o) \
               $(LIBRARY)(rotes2.o) \
               $(LIBRARY)(rtmi.o) \
               $(LIBRARY)(selfab.o) \
               $(LIBRARY)(setcur.o) \
               $(LIBRARY)(settrk.o) \
               $(LIBRARY)(sginit.o) \
               $(LIBRARY)(steep.o) \
               $(LIBRARY)(steeq.o) \
               $(LIBRARY)(steer.o) \
               $(LIBRARY)(stees.o) \
               $(LIBRARY)(stpair.o) \
               $(LIBRARY)(sub.o) \
               $(LIBRARY)(tdelay.o) \
               $(LIBRARY)(trac.o) \
               $(LIBRARY)(twob.o) \
               $(LIBRARY)(twoclu.o) 
