# Examen Parcial 1 - Computacion Paralela y Distribuida (CC3069)

CC      = gcc
CFLAGS  = -O2 -Wall -Wextra
OMPFLAGS= -fopenmp
BIN     = bin

PARALELOS   = $(BIN)/matrices $(BIN)/blur
SECUENCIALES= $(BIN)/matrices_original $(BIN)/blur_original

all: $(PARALELOS) $(SECUENCIALES)

$(BIN):
	mkdir -p $(BIN)

$(BIN)/matrices: paralelo/matrices_omp.c | $(BIN)
	$(CC) $(CFLAGS) $(OMPFLAGS) -o $@ $<

$(BIN)/blur: paralelo/blur_omp.c | $(BIN)
	$(CC) $(CFLAGS) $(OMPFLAGS) -o $@ $<

# Algoritmos base, tal como los trajo cada integrante.
$(BIN)/matrices_original: secuencial/matrices_secuencial.c | $(BIN)
	$(CC) $(CFLAGS) -o $@ $<

$(BIN)/blur_original: secuencial/blur_secuencial.c | $(BIN)
	$(CC) $(CFLAGS) -o $@ $<

bench: all
	./scripts/bench.sh

graficas:
	python3 scripts/graficar.py

clean:
	rm -rf $(BIN)

.PHONY: all bench graficas clean
