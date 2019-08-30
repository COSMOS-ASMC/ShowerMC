include  $(EPICSTOP)/site.config


objs = testKnock.o  epGetEffZA.o epReadMTbl.o epX0.o epCoulombC.o epStern.o  epExpot.o epKnock.o


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




