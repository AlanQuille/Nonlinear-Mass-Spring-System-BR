


CC      = g++
CXXLAGS = -Wall -std=gnu++11 -w

#add -O3 -Wa,-mbig-obj if using BDCSVD eigen class
all:: main

main: main.cpp
	$(CC) main.cpp $(CXXLAGS) -o main.o

clean:
	rm *.o
