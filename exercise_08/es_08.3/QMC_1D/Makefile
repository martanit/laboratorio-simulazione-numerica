CC = g++
CFLAGS = -O3 --std=c++11

qmc1d.exe : qmc1d.o random/random.o
	$(CC) random/random.o qmc1d.o -o qmc1d.exe
qmc1d.o : qmc1d.cpp
	$(CC) -c qmc1d.cpp -o qmc1d.o $(CFLAGS)
random.o : random/random.cpp random/random.h
	$(CC) -c random/random.cpp -o random/random.o $(CFLAGS)
clean :
	rm *.o ../../data/probability* random/*.o
run :
	./qmc1d.exe input.dat
