include  $(EPICSTOP)/site.config


objs = epBPfunc1.o  epBPAux.o  epwtSmpTbl.o getTotalXSofSeltzerBremsFunc.o epBrSeltzer.o epBPfuncH.o  epPairLowE.o epBrgeneric.o  epNormLPM.o epBPgeneini.o  epBPnormXs.o epBrgeneInte.o


getTotalXSofSeltzerBremsFunc.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




