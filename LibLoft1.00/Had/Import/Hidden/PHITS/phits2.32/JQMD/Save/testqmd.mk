include $(COSMOSTOP)/site.config
#all:
#	ifort -c -cpp  -I$COSMOSTOP/cosmos -I../Com/  qmd*f testqmd.f  mdp-uni.f  ../Com/{utl01.f,utl02.f,jamcomp.f,ndata01.f}  ../cjqmdout.f


#objs =   qmd*o testqmd.o  mdp-uni.o  ../Com/{utl01.o,utl02.o,jamcomp.o,ndata01.o}  ../cjqmdout.o
objs =   testqmd.o   ../Com/utl01.o ../Com/utl02.o ../Com/jamcomp.o ../Com/ndata01.o ../cjqmdout.o jqmd.o qmddflt.o

test$(ARCH): *.o
	$(LD)  -o $@ *.o -L$(DEST) -l$(LIBNAME) $(LDFLAGS)



