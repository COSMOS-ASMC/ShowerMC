include  $(EPICSTOP)/site.config


objs = testBremHx.o  epReadMTbl.o  psiim.o fbrem2.o epGetEffZA.o epX0Old.o epX0.o epCoulombC.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




