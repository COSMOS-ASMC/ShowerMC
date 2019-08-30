include  $(EPICSTOP)/site.config


objs =  epBPAux.o   DrawBrems.o epBrSeltzer.o epBPfuncH.o  epPairLowE.o epBrgeneric.o  epNormLPM.o epBPgeneini.o  

Drawbrems.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 


clean:;		@rm -f $(OBJS) core *~ a.out
