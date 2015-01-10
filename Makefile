.PHONY: all clean

CC=clang
CFLAGS=-lm

BIN = ising

all:$(BIN)

%.o: %.c
	$(CC) $(CFLAGS) -o $<
	
clean:
	rm -f $(BIN)
	rm -f *.o
	rm -f *~
