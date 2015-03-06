CC=c++
CFLAGS=-O3 -pthread -std=c++11

all: splice_detection

splice_detection: main.cpp
	$(CC) $(CFLAGS) main.cpp -o splice_detection
	#$(CC) main.o -o splice_detection

#main.o: main.cpp
#	$(CC) $(CFLAGS) main.cpp
clean:
	rm -rf splice_detection *.o
