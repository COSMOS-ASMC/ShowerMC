include $(COSMOSTOP)/site.config
#
#
vpath %.h  $(COSMOSINC)

vpath c%.f  $(COSMOSTOP)/Tracking/Test


src = testcdecayLeng.f cdecayLeng.f cdecayWEL.f

objs := $(src:.f=.o)


a.out:  $(objs)  
	$(LD) $(LDFLAGS) -o $@  $(objs) -L$(LIBLOC) -l$(LIBNAME)
###


