include $(LIBLOFT)/site.config

objs = cgetXsRatio.o cgetXsInterface2.o

a.out: $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)

