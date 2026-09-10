#include <stdio.h>

// Algoritmo 3: Multiplicación de matrices densas (C = A x B)
// Versión secuencial

// Tamaño de las matrices: 3 filas y 3 columnas
# define N 3


int main(void) {
    // Matriz A: 
    int matrizA[3][3] = {{5, 2, 1}, {2, 1, 2}, {4, 1, 3}};

    // Matriz B:
    int matrizB[3][3] = {{1, 4, 2}, {0, 3, 0}, {2, 1, 3}};

    // Matriz resultante C: (Donde se guardará el resultado de la multiplicación)
    int resultado[3][3];

    int n = 3;

    // Recorre las filas de la matriz A
    for (int i = 0; i < n; i++) {

        // Variable para almacenar la suma de los productos
        int suma = 0;

        // Recorre las columnas de la matriz B
        for (int j = 0; j < n; j++) {

            // Reinicia la variable suma para cada elemento 
            // de la matriz resultante
            suma = 0;

            // Recorre la fila i de A y la columna j de B
            for (int k = 0; k < n; k++) {

                // Multiplica los elementos correspondientes y acumula
                // los productos
                suma += matrizA[i][k] * matrizB[k][j];
            }
            // Guarda la suma obtenida en la posicion correspondiente
            resultado[i][j] = suma;
        }
    }

    // Imprime la matriz resultante C
    printf("Matriz resultado (C = A x B):\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d\t", resultado[i][j]);
        }
        printf("\n");
    }

    return 0;
}
