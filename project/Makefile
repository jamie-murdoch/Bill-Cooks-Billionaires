TARGET = edge_deletion
INCLUDES = -I C:/Program\ Files/IBM/ILOG/CPLEX_Studio1261/cplex/include/ilcplex/
LIBS = C:/Program\ Files/IBM/ILOG/CPLEX_Studio1261/cplex/bin/x64_win64/cplex1261.dll
CC = g++
CFLAGS =  -g -Wall -pedantic -Wno-long-long -ggdb -O3
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
