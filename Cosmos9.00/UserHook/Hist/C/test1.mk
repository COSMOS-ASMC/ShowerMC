objs = test1.o k90whist1.o

test1: $(objs)
	cc -o $@  $(objs) -lm  -lgsl -lgslcblas
