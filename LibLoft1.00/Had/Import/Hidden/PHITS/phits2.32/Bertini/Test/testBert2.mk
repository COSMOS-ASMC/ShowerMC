include $(LIBLOFT)/site.config
all:
	ifort -c -cpp  -I$(LIBLOFTINC) -I../../Com/  -Vaxlib -O0 -check all -traceback  -align dcommons ../cbert.f  ../../Com/ndata01.f   ../../Com/phitsrn.f ../../cphitsOut.f ../../cprePhits.f  ../../bert-bl0.f ../../bert-bl1.f ../../bert-bl2.f ../../ndata01blk.f ../../qmddflt.f

test$(ARCH):testBert.o
	$(LD)  -o $@ *.o -L$(DEST) -l$(LIBNAME) $(LDFLAGS)
