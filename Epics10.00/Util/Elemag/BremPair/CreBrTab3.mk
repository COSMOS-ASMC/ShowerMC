include  $(EPICSTOP)/site.config


objs = epBrSeltzer.o epCreBrTab.o epwtSmpTbl.o  epBPfuncH.o  epCreBrSTbH.o   epBPAux.o epCreBrSTbl1.o  epCrePrSTblH.o epCrePrSTbl1.o epBPfunc3.o  epPairLowE.o  epBrgeneric.o  epNormLPM.o epBPgeneini.o epPrgeneric.o epCreBrSTblS.o 
brem3.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:	@rm -f $(OBJS) core *~ a.out

include  $(EPICSTOP)/site.config






















drawbrems1.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out










