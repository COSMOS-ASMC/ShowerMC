include $(LIBLOFT)/site.config
#
#


src =  qgs01init.f  qgsjet01c.f testcol.f

objs := $(src:.f=.o)


testcol$(ARCH) :  $(objs)  
	$(LD) $(LDFLAGS)   -o  $@ $(objs) -L$(LIBLOC) -l$(LIBNAME)

