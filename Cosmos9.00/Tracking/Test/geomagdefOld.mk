include $(COSMOSTOP)/site.config
#

objs = testgeomagdef.o cmagneticDefOld.o


a.out:  $(objs)  
	$(LD) $(LDFLAGS) -o $@  $(objs) -L$(LIBLOC) -l$(LIBNAME)
