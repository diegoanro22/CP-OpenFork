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

## 2. Estrategia de paralelización

### 2.1 Matrices: el problema real no era la paralelización

El enunciado pregunta cómo asignar el trabajo *minimizando la cantidad de veces que los trabajadores
leen los mismos datos*. Al analizarlo encontramos que el cuello de botella del algoritmo original no
está en el reparto, sino en el orden de los ciclos.

En el orden `i-j-k` del original, el ciclo interno recorre `A[i][k]` avanzando de a un entero
(bien) pero `B[k][j]` avanzando de a `N` enteros (mal). Cada acceso a `B` cae en una línea de caché
distinta: el procesador trae 64 bytes de RAM y usa 4. Con N=1024 eso significa recorrer los 4 MB de
`B` completos por cada fila de `C` que se calcula.

Implementamos cuatro variantes para poder separar qué ganancia vino de dónde:

| Variante | Directiva | Idea |
|---|---|---|
| `seq` | ninguna | El triple ciclo original. Línea base. |
| `ijk` | `parallel for collapse(2)` sobre `i,j` | Paralelización ingenua: repartir las celdas de `C` y ya. |
| `ikj` | `parallel for schedule(static)` sobre `i` | Reordenar los ciclos a `i-k-j` para que todos los accesos sean contiguos. |
| `tiled` | `parallel for collapse(2)` sobre los bloques | `ikj` partido en bloques de 64×64 para reusar cada bloque desde caché. |

En `ikj`, el valor `A[i][k]` se carga una vez a un registro y se reusa durante todo el ciclo interno,
mientras `B[k][*]` y `C[i][*]` se barren de corrido. El patrón de acceso pasa a ser secuencial en las
tres matrices.

**Por qué se paraleliza `i` y nunca `k`.** En `ikj` el ciclo `k` acumula sobre `C[i][j]`: dos hilos
con distinto `k` escribirían la misma celda al mismo tiempo. Es una condición de carrera de libro y
para arreglarla habría que poner `atomic` en el ciclo más interno o una `reduction` sobre toda la
matriz `C`, y cualquiera de las dos costaría más de lo que ahorra. Paralelizando `i`, en cambio, cada
hilo es dueño exclusivo de un conjunto de filas de `C` y no hay nada que sincronizar.

En `tiled` pasa lo mismo un nivel más arriba: el `collapse(2)` va sobre los índices de bloque
`(ii, jj)` y no sobre `kk`, de modo que cada hilo se queda con un bloque completo de `C` y la
acumulación en `kk` ocurre secuencialmente dentro de ese bloque.

Sobre el **scheduling**: aquí el costo de cada iteración es idéntico —toda fila de `C` cuesta lo
mismo— así que `schedule(static)` es la opción correcta. Un `dynamic` solo agregaría el costo de
repartir trabajo en tiempo de ejecución sin ningún desbalance que corregir.

### 2.2 Blur: cómo se cortan los bordes de cada recorte

El enunciado pregunta cómo repartir la imagen *asegurándose de que se puedan calcular correctamente
los píxeles que quedan exactamente en los bordes de los recortes*. La respuesta corta es que en
memoria compartida ese problema **no existe**, y vale la pena explicar por qué.

Cortamos la imagen en **franjas horizontales de filas contiguas**, una por hilo, que es literalmente
lo que hace `#pragma omp parallel for schedule(static)` sobre el ciclo de filas: con 32 hilos y 7680
filas, cada hilo recibe un bloque de 240 filas seguidas.

Un hilo que procesa la última fila de su franja necesita leer la primera fila de la franja del
vecino. Esto sería un problema serio en memoria distribuida (habría que hacer intercambio de halos
por mensajes) y también lo sería si el filtro trabajara *in situ*, porque el vecino podría haber
sobrescrito esa fila antes de que la leamos. Pero como se **lee siempre de `A` y se escribe siempre
en `B`**, `A` nunca cambia durante el cálculo: cualquier hilo puede leer cualquier fila de `A` en
cualquier momento y siempre obtiene el valor original. El único borde que sí requiere tratamiento
especial es el de la imagen misma, y eso ya lo resolvía el recorte de ventana del secuencial.

La carrera de datos hay que buscarla en otro lado. En el secuencial original, las variables
`i0, i1, j0, j1, suma` y `contador` están declaradas al inicio de `main`, fuera de los ciclos
(líneas 5–6 de `blur_secuencial.c`). Si ese código se copia tal cual dentro de una región paralela,
OpenMP las trata como **compartidas** y los hilos se pisan el acumulador y los límites de la ventana
entre sí: el resultado sale mal, distinto en cada corrida y sin que el programa falle. Es el error
más fácil de cometer en este algoritmo. En la versión paralela el cálculo de un píxel se movió a la
función `pixel_promedio()`, donde todas esas variables son locales y por lo tanto privadas por
construcción.

Se probaron tres formas de repartir para poder justificar la elección con datos y no con intuición:
`filas` (`schedule(static)`), `collapse` (`collapse(2)`, que aplana los dos ciclos en uno de 33
millones de iteraciones) y `dynamic` (`schedule(dynamic, 64)`).

---
