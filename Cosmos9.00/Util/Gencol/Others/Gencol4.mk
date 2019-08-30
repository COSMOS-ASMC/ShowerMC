include  $(EPICSTOP)/site.config

objs =  Gencol4.o getDiffCode.o cmydecay.o

Gencol4$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos
clean:;		@rm -f $(OBJS) core *~ a.out

