include  $(EPICSTOP)/site.config

objs =  Gencol4.o  sibyllinit.o  sibyll2.1.o getDiffCode.o cmydecay.o

Gencol2$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos

clean:;		@rm -f $(OBJS) core *~ a.out

