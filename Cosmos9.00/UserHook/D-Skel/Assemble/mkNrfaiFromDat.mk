include $(COSMOSTOP)/site.config

objs = mkNrfaiFromDat.o

mkNrfai$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

