include $(LIBLOFT)/site.config

objs = testBert.o

cosmos$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)   $(LDFLAGS) 

