include  $(EPICSTOP)/site.config

objs = testScat.o epmulScat.o epReadMTbl.o  epGetEffZA.o epExpot.o epX0Old.o epX0.o epCoulombC.o epCompScrBr.o kexpi.o

a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 


clean:;		@rm -f $(OBJS) core *~ a.out







