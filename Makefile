CFLAGS=-O3

ngons: ngons.c
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f ngons
