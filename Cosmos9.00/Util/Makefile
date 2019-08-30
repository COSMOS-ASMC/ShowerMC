#
.PHONY: all clean drawGeomag showParam testPrim sampPrim drawRigCut


all: 
	@echo "---Use one of the following---" 
	@echo "make drawGeomag      For drawing Geomagnetism"
	@echo "make drawRigCut      For drawing rididity cutoff"
	@echo "make showParam       For making  showParam  "
	@echo "make testPrim        For testing primary data file"
	@echo "make sampPrim        For testing primary sampling"

drawGeomag:
	$(MAKE) -f drawGeomag.mk
drawRigCut:
	$(MAKE) -f drawRigCut.mk
showParam:
	$(MAKE) -f showParam.mk

testPrim:  
	$(MAKE) -f testPrim.mk
sampPrim:
	$(MAKE) -f sampPrim.mk

clean:
	@/bin/rm  -f *.o *.out  *.out.trace  core  *~ temp*  
###


