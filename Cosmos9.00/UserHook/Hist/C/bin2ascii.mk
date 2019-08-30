objs = bin2ascii.o k90whist1.o



bin2ascii: $(objs)
	cc -o $@  $(objs) -lm
