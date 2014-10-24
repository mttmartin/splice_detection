CC=cc
CFLAGS=-O3 -pthread

all: splice_detection

splice_detection: main.c
	$(CC) $(CFLAGS) main.c -o splice_detection
	#$(CC) main.o -o splice_detection

#main.o: main.c
#	$(CC) $(CFLAGS) main.c
clean:
	rm -rf splice_detection *.o
