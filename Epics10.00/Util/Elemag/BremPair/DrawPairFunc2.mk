include  $(EPICSTOP)/site.config

objs = epBPfunc2.o  epBPAux.o  epBPfuncH.o  epBPgeneini.o epPrgeneric.o  epPrgeneInte.o epBPnormXs.o   epNormLPM.o epBrSeltzer.o epwtSmpTbl.o  epPairLowE.o  DrawPair.o 



drawpair2.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out











