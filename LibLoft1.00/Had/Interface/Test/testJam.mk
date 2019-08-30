include $(LIBLOFT)/site.config

objs = testJam.o

testJam$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME) $(LDFLAGS)
