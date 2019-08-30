include  $(EPICSTOP)/site.config

objs = epmudsdvdr.o DrawdsdvdrAtv.o  epmuBPNgeneI.o epmuSetCnst.o epmuvmax.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




