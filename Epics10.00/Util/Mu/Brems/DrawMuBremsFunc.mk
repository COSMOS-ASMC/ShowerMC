include  $(EPICSTOP)/site.config



objs = epmuBrem.o  epmuAuxB.o  epwtSmpTbl.o DrawMuBrems.o  epMuBrEcheck.o



drawmubrems.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




