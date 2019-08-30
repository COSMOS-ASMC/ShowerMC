include $(LIBLOFT)/site.config

objs = test.o

test$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME) $(LDFLAGS)
