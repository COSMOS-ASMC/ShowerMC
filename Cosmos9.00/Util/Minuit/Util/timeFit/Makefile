include $(COSMOSTOP)/site.config
.PHONY:  all clean veryclean
all:
	make -f timeFit.mk
clean:
	@rm -f *.o core  *~ temp*.f a.out
veryclean: 
	@rm -f *.o core  *~  a.out temp*.f cosmos$(ARCH)
