include $(LIBLOFT)/site.config

objs = testcos2kf.o

testcos2kf$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME) $(LDFLAGS)
