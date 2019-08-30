include  $(EPICSTOP)/site.config


objs = epBPfunc3.o  epBPAux.o  epBPfuncH.o  epBPgeneini.o  epBrgeneric.o epBrgeneInte.o epBPnormXs.o  epNormLPM.o  epBrSeltzer.o    epPairLowE.o  epwtSmpTbl.o InteBrems.o 

Intebrems3.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




