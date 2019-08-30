include $(COSMOSTOP)/site.config

src = mkglauber.f

objs := $(src:.f=.o)


mkglauber$(ARCH):  $(objs) 
	$(LD) $(LDFLAGS) -o $@  $(objs) -L$(LIBLOC) -l$(LIBNAME)  -L$(LIBLOFT)/lib/$(ARCH) -lloft

clean:
	rm -f a.out *.o
