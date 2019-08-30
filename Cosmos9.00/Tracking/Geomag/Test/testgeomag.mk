include $(COSMOSTOP)/site.config

objs = testgeomag.o 

testgeomag$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)


