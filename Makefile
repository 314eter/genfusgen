CFLAGS=-O3

genfusgen: genfusgen.c
	$(CC) $(CFLAGS) $^ -o $@

.PHONY: clean

clean:
	rm -f genfusgen
