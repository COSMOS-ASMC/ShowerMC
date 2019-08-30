include $(LIBLOFT)/site.config

objs = xsecFit.o

execFit: $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs)  -L$(DEST) -l$(LIBNAME) -lminuit 
