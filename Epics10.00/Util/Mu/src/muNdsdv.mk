include  $(EPICSTOP)/site.config

objs = epmuNdsdv.o DrawMuNdsdv.o  epmuBPNgeneI.o epmuSetCnst.o  epmuAuxNBB.o 

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




	