include  $(EPICSTOP)/site.config


objs =  DrawMuNdsdv.o epMuNiEcheck.o epmuAuxNBB.o epmuNdsdv.o


drawmunuci.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




