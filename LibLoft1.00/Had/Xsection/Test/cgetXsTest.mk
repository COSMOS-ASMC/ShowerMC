include $(LIBLOFT)/site.config

objs = cgetXsTest.o cgetXsInterface2.o

a.out: $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)

