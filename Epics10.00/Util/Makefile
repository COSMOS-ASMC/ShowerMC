.PHONY: all drawConfig test


all: 
	@echo "To compile drawConfig.f;  make drawConfig"
drawConfig:
	make -f drawConfig.mk
	rm -f /tmp/*.f
test: 
	make -f test.mk
	rm -f /tmp/*.f

clean:
	@rm -f *.o core  *~  temp*.f

veryclean:
	@rm -f *.o core  *~  temp*.f*  drawconfig

rmtemp: 
	@rm  -f temp*.f

