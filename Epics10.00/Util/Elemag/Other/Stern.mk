include  $(EPICSTOP)/site.config

objs = epReadMTbl.o testStern.o  epExpot.o epStern.o epwtStern.o epGetEffZA.o epX0.o epCoulombC.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




