include  $(EPICSTOP)/site.config

objs =  main.o  readMedia.o setBlockData.o init.o

showXsec$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -v   -L$(DEST) -l$(LIBNAME) -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos 

