TARGET = subtour
INCLUDES = -I /Users/lawson/Applications/IBM/ILOG/CPLEX_Studio1261/cplex/include/ilcplex/
LIBS = /Users/lawson/Applications/IBM/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_osx/static_pic/libcplex.a
CC = g++
CFLAGS =  -g -Wall -pedantic -Wno-long-long -ggdb 
#-O3
#-O3
.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)