include $(COSMOSTOP)/site.config

objs = translate.o catmncTcos.o

translate$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS) -L$(LIBLOFT)/lib/$(ARCH) -lloft
