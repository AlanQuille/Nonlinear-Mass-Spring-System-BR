
CC      = g++
CXXLAGS = -Wall -std=gnu++11 -w
CFLAG = -m64

#add -O3 -Wa,-mbig-obj if using BDCSVD eigen class
all:: main

main: main.cpp
	$(CC) main.cpp $(CXXLAGS) -I"C:/Program Files/MATLAB/R2017b/bin/win64" -L"C:/Program Files/MATLAB/R2017b/extern/lib/win64/mingw64" -llibmat -llibmx -o main.o

clean:
	rm *.o
