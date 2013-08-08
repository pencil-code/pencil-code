all: 1 2 3 4 5

%: cycle.h sz.h main.c %.c
	gcc main.c $@.c -o $@ -Wall -O3

cycle.h:
	./mkcycle

clean:
	-rm *~ cycle.h ? out
