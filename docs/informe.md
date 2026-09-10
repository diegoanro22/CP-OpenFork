# Informe técnico

**CC3069 — Computación Paralela y Distribuida · Examen Parcial 1**
Consultora OpenFork U · Diego Rosales, Jose Lopez, Oliver Viau

---

## 1. Contexto y datos

De los tres problemas que traía el equipo elegimos dos: la **multiplicación de matrices densas**
(problema 3) y el **filtro de desenfoque sobre imagen satelital** (problema 5). El tercero quedó
fuera porque los dos elegidos plantean preguntas de optimización distintas entre sí —uno es un
problema de reuso de datos en caché y el otro de repartición de un dominio 2D con dependencias de
vecindad— y eso daba más material que dos problemas parecidos.

### 1.1 Multiplicación de matrices

La propuesta secuencial original está en [`secuencial/matrices_secuencial.c`](../secuencial/matrices_secuencial.c):
un triple ciclo `i-j-k` sobre matrices de 3×3 escritas a mano en el código. Calcula bien, pero para
medir cualquier cosa hace falta un tamaño donde el trabajo se note.

En la versión de [`paralelo/matrices_omp.c`](../paralelo/matrices_omp.c) las matrices son cuadradas
de lado `N` configurable y se llenan con `rand()` bajo una semilla fija (`20260907`). La semilla fija
es lo que permite que los tres integrantes midan sobre exactamente los mismos datos en máquinas
distintas y que los tiempos sean comparables.

Dos decisiones sobre la estructura en memoria:

- **Enteros, no punto flotante.** Los valores van de 0 a 9, así cada celda del resultado se queda muy
  por debajo del límite de `int` incluso con N=4096. La razón de fondo es la verificación: como la
  suma de enteros sí es asociativa, el checksum del resultado tiene que salir **idéntico** en las
  cuatro variantes. Si hubiéramos usado `double`, un checksum distinto podría deberse tanto a una
  condición de carrera como al simple reordenamiento de la suma, y no habría forma de distinguirlas.
- **Una sola reserva contigua** (`int *` indexado como `A[i*N+j]`) en lugar del arreglo de punteros
  del original. Con `malloc` por fila, las filas quedan dispersas en el heap y el prefetcher del
  procesador no puede anticipar nada; con un bloque contiguo, recorrer una fila es recorrer memoria
  consecutiva. Esta decisión sola ya cambia el rendimiento antes de tocar OpenMP.

El tamaño de referencia para las mediciones es **N=1024**, con una corrida adicional en N=2048 para
verificar el comportamiento a mayor escala.

### 1.2 Filtro de desenfoque

El secuencial original está en [`secuencial/blur_secuencial.c`](../secuencial/blur_secuencial.c).
Lee las dimensiones y los píxeles desde entrada estándar y, para cada píxel, promedia su valor con
el de sus vecinos en una ventana de 3×3 recortada contra los bordes de la imagen.

La versión paralela ([`paralelo/blur_omp.c`](../paralelo/blur_omp.c)) conserva el kernel tal cual y
cambia únicamente de dónde salen los datos: la imagen se genera en memoria con la misma semilla fija,
con el tamaño 8K del enunciado (**7680 × 4320 = 33.2 millones de píxeles**, 33 MB por buffer).

Se generan los píxeles en vez de leer un archivo por dos razones. La primera es que un PGM de 8K pesa
33 MB y no tiene sentido meterlo en el repositorio. La segunda es que leer 33 MB de disco toma más
tiempo que el filtro completo, así que el I/O ahogaría justamente lo que queremos medir.

Los píxeles son `unsigned char` en un bloque contiguo, y hay **dos buffers separados**: se lee de `A`
y se escribe en `B`. Esto no es un detalle de implementación, es la decisión que hace que todo el
algoritmo sea paralelizable sin sincronización, y se explica en la sección 2.2.

---
