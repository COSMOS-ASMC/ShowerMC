include $(COSMOSTOP)/site.config

.PHONY: all clean initObs



all:
	@echo "--- Use one of following---:"
	@echo "make initObs         For testing initObs"

initObs:  
	$(MAKE) -f initObs.mk initObs

clean:
	/bin/rm  -f *.o *.out  *.out.trace  core  *~ temp*

