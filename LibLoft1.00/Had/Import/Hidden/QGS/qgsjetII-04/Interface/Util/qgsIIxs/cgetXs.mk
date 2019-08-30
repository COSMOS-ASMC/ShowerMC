include $(LIBLOFT)/site.config

objs = cgetXs.o cgetXsInterface.o

a.out: $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)

