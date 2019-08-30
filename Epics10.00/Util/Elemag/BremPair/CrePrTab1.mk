include  $(EPICSTOP)/site.config

objs = epBrSeltzer.o epCrePrTab.o epwtSmpTbl.o epBPfuncH.o  epCreBrSTbH.o  epBPAux.o epCreBrSTbl1.o  epCrePrSTblH.o epCrePrSTbl1.o epBPfunc1.o epPairLowE.o epPrgeneric.o epNormLPM.o epBPgeneini.o epBrgeneric.o epPrgeneInte.o  epBrgeneInte.o epBPnormXs.o


pair1.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out







