include  $(EPICSTOP)/site.config

objs =  DrawMuNdEdx.o  epmuBPNgeneI.o epmuAuxNBB.o epmuSetCnst.o epmuNdsdv.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




