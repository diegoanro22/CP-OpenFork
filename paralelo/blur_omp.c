// Algoritmo 5: filtro de desenfoque sobre una imagen, con OpenMP.
// Uso: ./blur <seq|filas|collapse|dynamic> <M> <N>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define SEMILLA 20260907

static unsigned char *reservar(int m, int n) {
    unsigned char *p = malloc((size_t)m * (size_t)n);
    if (p == NULL) {
        fprintf(stderr, "Error de memoria (%d x %d).\n", m, n);
        exit(1);
    }
    return p;
}

static void generar(unsigned char *a, int m, int n) {
    for (size_t i = 0; i < (size_t)m * (size_t)n; i++) {
        a[i] = (unsigned char)(rand() % 256);
    }
}

static long long checksum(const unsigned char *b, int m, int n) {
    long long s = 0;
    for (size_t i = 0; i < (size_t)m * (size_t)n; i++) {
        s += b[i];
    }
    return s;
}

// Mismo calculo que el secuencial. Las temporales van locales para que sean
// privadas; en el original estaban declaradas afuera, en main.
static inline unsigned char pixel_promedio(const unsigned char *A, int m, int n,
                                           int i, int j) {
    int i0 = (i - 1 > 0) ? i - 1 : 0;
    int i1 = (i + 1 < m - 1) ? i + 1 : m - 1;
    int j0 = (j - 1 > 0) ? j - 1 : 0;
    int j1 = (j + 1 < n - 1) ? j + 1 : n - 1;

    int suma = 0;
    int contador = 0;

    for (int r = i0; r <= i1; r++) {
        for (int c = j0; c <= j1; c++) {
            suma += A[(size_t)r * n + c];
            contador++;
        }
    }
    return (unsigned char)(suma / contador);
}

static void blur_seq(const unsigned char *A, unsigned char *B, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            B[(size_t)i * n + j] = pixel_promedio(A, m, n, i, j);
        }
    }
}

// schedule(static) le da a cada hilo una franja de filas contiguas.
// Los bordes de cada franja no necesitan nada extra: se lee de A y se escribe
// en B, asi que A nunca cambia mientras se calcula.
static void blur_filas(const unsigned char *A, unsigned char *B, int m, int n) {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            B[(size_t)i * n + j] = pixel_promedio(A, m, n, i, j);
        }
    }
}

// Aplana los dos ciclos en uno solo de m*n iteraciones.
static void blur_collapse(const unsigned char *A, unsigned char *B, int m, int n) {
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            B[(size_t)i * n + j] = pixel_promedio(A, m, n, i, j);
        }
    }
}

// Reparto por demanda, para comparar contra el static.
static void blur_dynamic(const unsigned char *A, unsigned char *B, int m, int n) {
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            B[(size_t)i * n + j] = pixel_promedio(A, m, n, i, j);
        }
    }
}

int main(int argc, char **argv) {
    const char *variante = (argc > 1) ? argv[1] : "seq";
    int m = (argc > 2) ? atoi(argv[2]) : 7680;   // 8K por defecto
    int n = (argc > 3) ? atoi(argv[3]) : 4320;

    if (m <= 0 || n <= 0) {
        fprintf(stderr, "Dimensiones invalidas.\n");
        return 1;
    }

    unsigned char *A = reservar(m, n);
    unsigned char *B = reservar(m, n);

    srand(SEMILLA);
    generar(A, m, n);

    double t0 = omp_get_wtime();

    if (strcmp(variante, "seq") == 0) {
        blur_seq(A, B, m, n);
    } else if (strcmp(variante, "filas") == 0) {
        blur_filas(A, B, m, n);
    } else if (strcmp(variante, "collapse") == 0) {
        blur_collapse(A, B, m, n);
    } else if (strcmp(variante, "dynamic") == 0) {
        blur_dynamic(A, B, m, n);
    } else {
        fprintf(stderr, "Variante desconocida: %s (use seq|filas|collapse|dynamic)\n", variante);
        return 1;
    }

    double t1 = omp_get_wtime();

    printf("algoritmo=blur variante=%s N=%dx%d hilos=%d tiempo=%.6f checksum=%lld\n",
           variante, m, n,
           (strcmp(variante, "seq") == 0) ? 1 : omp_get_max_threads(),
           t1 - t0, checksum(B, m, n));

    free(A);
    free(B);
    return 0;
}
