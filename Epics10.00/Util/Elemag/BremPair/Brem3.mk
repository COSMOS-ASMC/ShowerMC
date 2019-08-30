include  $(EPICSTOP)/site.config

objs = epBPfunc3.o testBrem.o epX0.o epCoulombC.o epX0Old.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




