include $(COSMOSTOP)/site.config
.PHONY:  all clean veryclean
all:

	c++ -c -w ../cmain.cc -I../../include
	c++ -c -w chookHybAS.cc -I../../include
	c++ -c -w ../ctemplCeren.cc -I../../include
	c++ -c -w chook.cc -I../../include
	 /lib/cpp -C -P -traditional -I$(COSMOSINC) ../cmymain.f > cmymain.f
	ifort -c -w cmymain.f
	c++ -o cosmos$(ARCH) cmain.o cmymain.o chook.o chookHybAS.o ctemplCeren.o -L$(COSMOSTOP)/lib/$(ARCH) -lcosmos `cat ifclib `
clean:
	@rm -f *.o core  *~ temp*.f a.out
veryclean: 
	@rm -f *.o core  *~  a.out temp*.f cosmos$(ARCH)

