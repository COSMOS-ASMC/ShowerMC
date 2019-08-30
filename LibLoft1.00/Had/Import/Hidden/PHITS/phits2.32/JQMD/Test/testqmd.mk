include $(LIBLOFT)/site.config
all:
#	ifort -c -cpp  -I$(LIBLOFT)/Headeras -I../Com/  -Vaxlib -check all -traceback  -align dcommons  jqmd.f testqmd.f  mdp-uni.f   ../Com/{jamcomp.f,ndata01.f,phitsrn.f}  ../cphitsOut.f ../cprePhits.f  ../cphits2cos.f
objs = testqmd.o


test$(ARCH): testqmd.o
	$(LD)  -o $@ *.o -L$(DEST) -l$(LIBNAME) $(LDFLAGS)
	rm testqmd.o
