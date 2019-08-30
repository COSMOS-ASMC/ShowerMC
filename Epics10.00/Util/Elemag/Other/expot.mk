include  $(EPICSTOP)/site.config

objs = epReadMTbl.o testExpot.o  epExpot.o epGetEffZA.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




