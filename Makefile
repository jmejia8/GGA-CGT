include Make.inc
CFLAGS-add += -DMATHLIB_STANDALONE

SRC = $(wildcard *.c)
TAR = $(SRC:.c=.so)

# all: $(TAR)

GGA-CGT: GGA-CGT.CPP
	$(CC) $(CFLAGS-add) -o  $@ $^

clean:
	rm GGA-CGT