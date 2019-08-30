include  $(EPICSTOP)/site.config


objs = testPhotoE.o  epGetEffZA.o epReadMTbl.o epX0.o epCoulombC.o epphotoE.o epSetPhotoE.o


a.out: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 



clean:;		@rm -f $(OBJS) core *~ a.out




