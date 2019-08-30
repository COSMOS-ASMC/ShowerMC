include $(COSMOSTOP)/site.config

.PHONY: all clean initObs



all:
	@echo "--- Use one of following---:"
	@echo "make mkInc         For testing cmkIncident"
	@echo "make sPrimAng      For testing csPrimAng"

mkInc:  
	$(MAKE) -f mkIncident.mk mkInc

sPrimAng:
	$(MAKE) -f sPrimAng.mk sPrimAng
clean:
	/bin/rm  -f *.o *.out  *.out.trace  core  *~ temp*
