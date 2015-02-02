.PHONY: all clean

CC=clang
CFLAGS=-lm -Wall -Wextra -O3

BIN = ising

all:$(BIN)

%.o: %.c
	$(CC) $(CFLAGS) -o $<
	
clean:
	rm -f $(BIN)
	rm -f *.o
	rm -f *~
