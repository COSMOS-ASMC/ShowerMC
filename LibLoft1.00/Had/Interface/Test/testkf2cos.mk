include $(LIBLOFT)/site.config

objs = testkf2cos.o

testkf2cos$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME) $(LDFLAGS)
