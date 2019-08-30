include $(LIBLOFT)/site.config

objs = test.o 


sofia$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)


test.o: $(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zptcl.h 

