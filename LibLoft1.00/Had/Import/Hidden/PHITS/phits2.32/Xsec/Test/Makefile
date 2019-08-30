include $(COSMOSTOP)/site.config
all:
	make -f testXsec.mk
	make -f testXsec.mk test$(ARCH)
clean:
	@rm -f $(OBJS) core *~

veryclean:
	@rm -f $(OBJS) core *~ temp*.f

