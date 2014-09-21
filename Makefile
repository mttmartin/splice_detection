CC=cc
CFLAGS=-O3 -c

all: splice_detection

splice_detection: main.o
	$(CC) main.o -o splice_detection

main.o: main.c
	$(CC) $(CFLAGS) main.c
clean:
	rm -rf splice_detection *.o
