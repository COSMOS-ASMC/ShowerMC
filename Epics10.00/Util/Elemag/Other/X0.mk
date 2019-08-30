include  $(EPICSTOP)/site.config

objs = epReadMTbl.o testX0.o  epX0.o epGetEffZA.o epCoulombC.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




