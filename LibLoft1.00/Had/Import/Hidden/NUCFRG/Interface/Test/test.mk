include $(LIBLOFT)/site.config

objs = testsampNucFrag.o

sampNucFrag$(ARCH): $(objs)
	$(LD)  -o $@ $(objs) -L$(DEST) -l$(LIBNAME)  $(LDFLAGS)


testsampNucFrag.o: $(LIBLOFTINC)/Zptcl.h 
