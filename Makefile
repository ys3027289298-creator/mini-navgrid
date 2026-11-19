CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: mini-navgrid

mini-navgrid: src/navgrid.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: mini-navgrid
	./mini-navgrid

clean:
	rm -f mini-navgrid
