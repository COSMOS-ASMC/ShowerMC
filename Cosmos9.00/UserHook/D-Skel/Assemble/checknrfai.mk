include $(COSMOSTOP)/site.config

objs = checknrfai.o

checknrfai$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)

