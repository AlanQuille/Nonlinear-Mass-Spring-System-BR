


CC      = g++
CXXLAGS = -Wall


all:: main

main: main.cpp
	$(CC) main.cpp $(CXXLAGS) -o main.o

clean:
	rm *.o

