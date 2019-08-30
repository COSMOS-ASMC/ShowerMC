#  to be used at tastim599
#  #define USECOSMOS
#   must be used at the top of test1.c
objs = test1.o k90whist1.o

test1: $(objs)
	cc -o $@  $(objs) -lm  -L/DB/TAMCDB/F/lib/PCLinuxIFC -lcosmos `cat /DB/TAMCDB/F/ifclib`
