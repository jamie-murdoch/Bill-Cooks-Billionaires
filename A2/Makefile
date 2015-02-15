SUBTOUROBJS= subtour_hw.o lp.o util.o
CPLEX_LIB = /software/cplex/distribution/cplex/lib/x86-64_linux/static_pic/libcplex.a
CPLEX_INCLUDE = /software/cplex/distribution/cplex/include/ilcplex

CC= gcc
CFLAGS=  -O3 -g -Wall -pedantic -Wno-long-long -I$(CPLEX_INCLUDE)


subtour: $(SUBTOUROBJS)
	  $(CC) $(CFLAGS) -o subtour $(SUBTOUROBJS) $(CPLEX_LIB) -lm -lpthread

clean:
	rm -f *.o subtour

subtour_hw.o: subtour_hw.c util.h
lp.o: lp.c lp.h
util.o: util.c util.h
