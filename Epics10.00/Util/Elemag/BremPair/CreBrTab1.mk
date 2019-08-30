include  $(EPICSTOP)/site.config


objs = epBPfunc1.o  epBPAux.o  epBPfuncH.o  epBPgeneini.o  epBrgeneric.o epBrgeneInte.o epBPnormXs.o  epNormLPM.o  epBrSeltzer.o    epPairLowE.o  epwtSmpTbl.o  epCreBrSTbH.o    epCreBrSTbl1.o    epCreBrSTblS.o  epCreBrTab.o 


brem1.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH)   -lcosmos


clean:;		@rm -f $(OBJS) core *~ a.out























