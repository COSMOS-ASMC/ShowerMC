include $(LIBLOFT)/site.config

objs = Bsin.o

Bsin$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)

