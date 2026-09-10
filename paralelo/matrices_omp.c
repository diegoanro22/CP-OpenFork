// Algoritmo 3: multiplicacion de matrices densas (C = A x B) con OpenMP.
// Uso: ./matrices <seq|ijk|ikj|tiled> <N>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

// Semilla fija para que todos midan sobre los mismos datos.
#define SEMILLA 20260907

#define BLOQUE 64

static int min_int(int a, int b) {
    return (a < b) ? a : b;
}

// Un solo bloque contiguo en vez de un arreglo de punteros por fila.
static int *reservar(int n) {
    int *m = malloc((size_t)n * (size_t)n * sizeof(int));
    if (m == NULL) {
        fprintf(stderr, "Error de memoria (N=%d).\n", n);
        exit(1);
    }
    return m;
}

// Valores de 0 a 9 para que el resultado no se pase de int.
static void generar(int *m, int n) {
    for (int i = 0; i < n * n; i++) {
        m[i] = rand() % 10;
    }
}

// Sirve para comparar que todas las variantes den lo mismo.
static long long checksum(const int *c, int n) {
    long long s = 0;
    for (int i = 0; i < n * n; i++) {
        s += c[i];
    }
    return s;
}

// Version original, sin directivas. Es la linea base.
static void mult_seq(const int *A, const int *B, int *C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int suma = 0;
            for (int k = 0; k < n; k++) {
                suma += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = suma;
        }
    }
}

// Reparte las celdas de C entre los hilos sin tocar el orden de los ciclos.
static void mult_ijk(const int *A, const int *B, int *C, int n) {
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // suma y k van adentro para que queden privadas.
            int suma = 0;
            for (int k = 0; k < n; k++) {
                suma += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = suma;
        }
    }
}

// Orden i-k-j: los tres accesos quedan contiguos en memoria.
// Se paraleliza i y no k, porque k acumula sobre C[i][j].
static void mult_ikj(const int *A, const int *B, int *C, int n) {
    memset(C, 0, (size_t)n * (size_t)n * sizeof(int));

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            int a = A[i * n + k];
            for (int j = 0; j < n; j++) {
                C[i * n + j] += a * B[k * n + j];
            }
        }
    }
}

// ikj por bloques, para reusar cada bloque desde cache.
// El collapse va sobre (ii, jj): cada hilo se queda con un bloque entero de C.
static void mult_tiled(const int *A, const int *B, int *C, int n) {
    memset(C, 0, (size_t)n * (size_t)n * sizeof(int));

    #pragma omp parallel for collapse(2) schedule(static)
    for (int ii = 0; ii < n; ii += BLOQUE) {
        for (int jj = 0; jj < n; jj += BLOQUE) {
            for (int kk = 0; kk < n; kk += BLOQUE) {

                int i_max = min_int(ii + BLOQUE, n);
                int j_max = min_int(jj + BLOQUE, n);
                int k_max = min_int(kk + BLOQUE, n);

                for (int i = ii; i < i_max; i++) {
                    for (int k = kk; k < k_max; k++) {
                        int a = A[i * n + k];
                        for (int j = jj; j < j_max; j++) {
                            C[i * n + j] += a * B[k * n + j];
                        }
                    }
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    const char *variante = (argc > 1) ? argv[1] : "seq";
    int n = (argc > 2) ? atoi(argv[2]) : 1024;

    if (n <= 0) {
        fprintf(stderr, "N invalido.\n");
        return 1;
    }

    int *A = reservar(n);
    int *B = reservar(n);
    int *C = reservar(n);

    srand(SEMILLA);
    generar(A, n);
    generar(B, n);

    double t0 = omp_get_wtime();

    if (strcmp(variante, "seq") == 0) {
        mult_seq(A, B, C, n);
    } else if (strcmp(variante, "ijk") == 0) {
        mult_ijk(A, B, C, n);
    } else if (strcmp(variante, "ikj") == 0) {
        mult_ikj(A, B, C, n);
    } else if (strcmp(variante, "tiled") == 0) {
        mult_tiled(A, B, C, n);
    } else {
        fprintf(stderr, "Variante desconocida: %s (use seq|ijk|ikj|tiled)\n", variante);
        return 1;
    }

    double t1 = omp_get_wtime();

    printf("algoritmo=matrices variante=%s N=%d hilos=%d tiempo=%.6f checksum=%lld\n",
           variante, n,
           (strcmp(variante, "seq") == 0) ? 1 : omp_get_max_threads(),
           t1 - t0, checksum(C, n));

    free(A);
    free(B);
    free(C);
    return 0;
}
