


CC      = g++
CXXLAGS = -Wall -O3 -Wa,-mbig-obj -std=gnu++11 -w


all:: main

main: main.cpp
	$(CC) main.cpp $(CXXLAGS) -o main.o

clean:
	rm *.o
