include $(COSMOSTOP)/site.config
#

#  objs = testgeomagdef.o cmagneticDefOld.o
objs = testgeomagdef.o 


a.out:  $(objs)  
	$(LD) $(LDFLAGS) -o $@  $(objs) -L$(LIBLOC) -l$(LIBNAME)
