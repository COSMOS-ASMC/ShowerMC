include $(LIBLOFT)/site.config
#
#


src =  testcol.f 

objs := $(src:.f=.o)


testcol$(ARCH) :  $(objs)  
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(LIBLOC) -l$(LIBNAME)

###


testcol.o: $(LIBLOFTINC)/Zptcl.h \
	$(LIBLOFTINC)/Zcode.h \
	$(LIBLOFTINC)/Zevhnv.h
