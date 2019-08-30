include  $(EPICSTOP)/site.config

objs = epmudsdvdr.o  epmuBPNgeneI.o epmuSetCnst.o epmuvmax.o epmuPrvVsdEdx.o epmuAuxP.o epmudsdv.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




