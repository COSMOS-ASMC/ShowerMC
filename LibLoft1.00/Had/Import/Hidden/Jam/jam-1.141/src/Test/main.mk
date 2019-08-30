include $(LIBLOFT)/site.config

objs = main.o mydummy.o

cosmos$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)


