#include <stdio.h>
#include <stdlib.h>

int main(void) {
    int M, N, i, j, r, c;
    int suma, contador, i0, i1, j0, j1;
    unsigned char **A; // matriz original
    unsigned char **B; // matriz resultado

    /* ---- Leer dimensiones de la imagen ---- */
    printf("Filas (M) y columnas (N) de la imagen: ");
    if (scanf("%d %d", &M, &N) != 2 || M <= 0 || N <= 0) {
        fprintf(stderr, "Dimensiones invalidas.\n");
        return 1;
    }

    /* ---- Reservar matrices A y B ---- */
    A = malloc((size_t)M * sizeof(unsigned char *));
    B = malloc((size_t)M * sizeof(unsigned char *));
    if (A == NULL || B == NULL) {
        fprintf(stderr, "Error de memoria.\n");
        return 1;
    }
    for (i = 0; i < M; i++) {
        A[i] = malloc((size_t)N * sizeof(unsigned char));
        B[i] = malloc((size_t)N * sizeof(unsigned char));
        if (A[i] == NULL || B[i] == NULL) {
            fprintf(stderr, "Error de memoria.\n");
            return 1;
        }
    }

    /* ---- Cargar los pixeles de A[M][N] (valores 0..255) ---- */
    printf("Ingrese los %d x %d valores de pixel (0-255):\n", M, N);
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            int valor;
            scanf("%d", &valor);
            A[i][j] = (unsigned char)valor;
        }
    }

    /* ---- Recorrer cada pixel y calcular su nuevo valor ---- */
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {

            /* Ventana valida para A[i][j], recortada a los bordes */
            i0 = (i - 1 > 0) ? i - 1 : 0;
            i1 = (i + 1 < M - 1) ? i + 1 : M - 1;
            j0 = (j - 1 > 0) ? j - 1 : 0;
            j1 = (j + 1 < N - 1) ? j + 1 : N - 1;

            suma = 0;
            contador = 0;

            for (r = i0; r <= i1; r++) {
                for (c = j0; c <= j1; c++) {
                    suma += A[r][c];
                    contador++;
                }
            }

            B[i][j] = (unsigned char)(suma / contador);
        }
    }

    /* ---- Guardar o mostrar B[M][N] ---- */
    printf("\nImagen resultante (B):\n");
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            printf("%3d ", B[i][j]);
        }
        printf("\n");
    }

    /* ---- Liberar memoria ---- */
    for (i = 0; i < M; i++) {
        free(A[i]);
        free(B[i]);
    }
    free(A);
    free(B);

    return 0;
}
