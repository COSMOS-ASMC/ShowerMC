include $(LIBLOFT)/site.config

objs = testnela.o

cosmos$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME) -lloft  $(LDFLAGS) 
