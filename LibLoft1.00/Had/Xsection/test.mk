include $(COSMOSTOP)/site.config

objs = cpdgXs.o

a.out: $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)


