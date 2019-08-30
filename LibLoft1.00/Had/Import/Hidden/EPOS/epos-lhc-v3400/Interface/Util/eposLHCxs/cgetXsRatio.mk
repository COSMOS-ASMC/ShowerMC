include $(LIBLOFT)/site.config

objs = cgetXsRatio.o cgetXsInterface.o

a.out: $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)

