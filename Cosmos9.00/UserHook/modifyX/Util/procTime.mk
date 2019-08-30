include $(COSMOSTOP)/site.config

#  objs = chook.o interface.o crbinPortion.o cxppipe.o kxplcy.o cminTime2WebSec.o
objs = procTime0.o procTime.o 

procTime$(ARCH): $(objs)
	$(LD) $(LDFLAGS) -o $@ $(objs) -L$(DEST) -l$(LIBNAME)


chook.o: ZtimeAna.h \
	ZprivateSub.h

