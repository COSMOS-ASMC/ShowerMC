include  $(EPICSTOP)/site.config


objs = epBPfunc1.o  epBPAux.o  epwtSmpTbl.o  epBrSeltzer.o epBPfuncH.o  epPairLowE.o epPrgeneric.o  epNormLPM.o epBPgeneini.o epPrgeneInte.o epBPnormXs.o IntePair.o



IntePair1.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




