include  $(EPICSTOP)/site.config

objs = epBPfunc1.o testBremComposit.o epX0.o epCoulombC.o epX0Old.o epBPAux.o epGetEffZA.o epReadMTbl.o


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




