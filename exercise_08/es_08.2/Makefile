CC = g++
CFLAGS = -O3 --std=c++11

main.exe : main.o random/random.o
	$(CC) random/random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random/random.cpp random/random.h
	$(CC) -c random/random.cpp -o random/random.o $(CFLAGS)
clean :
	rm *.o ../data/* random/*.o
run :
	./main.exe
