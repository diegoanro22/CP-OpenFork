# OpenFork U

Examen Parcial 1 de Computación Paralela y Distribuida (CC3069), UVG, Semestre II 2026.

**Integrantes:** Diego Rosales · Jose Lopez · Oliver Viau

De los tres problemas que traía el equipo elegimos dos para paralelizar con OpenMP:

- **Multiplicación de matrices densas** (problema 3), porque el cuello de botella no está en repartir
  el trabajo sino en cómo se lee la memoria, y eso permitía atacar el reuso de datos en caché.
- **Filtro de desenfoque sobre imagen satelital** (problema 5), porque plantea el problema de cortar
  un dominio 2D entre trabajadores cuando cada píxel depende de sus vecinos.

El informe completo, con la estrategia de paralelización y las mediciones de speedup y eficiencia,
está en [`docs/informe.md`](docs/informe.md).

## Estructura

```
secuencial/   Los algoritmos base, tal como los trajo cada integrante.
paralelo/     Las versiones con OpenMP, con varias estrategias cada una.
scripts/      Benchmark de speedup y eficiencia.
docs/         Informe, mediciones y evidencia de las corridas.
```

## Cómo correrlo

Hace falta gcc con soporte de OpenMP.

```sh
make                        # compila todo en bin/
./bin/matrices ikj 1024     # variantes: seq | ijk | ikj | tiled
./bin/blur filas 7680 4320  # variantes: seq | filas | collapse | dynamic
```

Cada binario imprime el tiempo y un checksum del resultado. El checksum tiene que ser el mismo en
todas las variantes del mismo algoritmo; si cambia, hay una condición de carrera.

Para reproducir las mediciones completas:

```sh
./scripts/bench.sh tu-nombre
```

Barre de 1 hasta la cantidad de hilos de la máquina, corre cada configuración tres veces y se queda
con la mediana, verifica los checksums contra la versión secuencial y deja los resultados en
`docs/resultados/`. Las gráficas de speedup y eficiencia salen de ahí:

```sh
python3 scripts/graficar.py tu-nombre   # requiere matplotlib
```

En una máquina lenta el benchmark se puede acortar:

```sh
REPS=1 HILOS="1 2 4" N_MAT=512 M_IMG=1920 N_IMG=1080 ./scripts/bench.sh tu-nombre
```

Los dos algoritmos originales también se pueden correr tal cual venían:

```sh
./bin/matrices_original
echo "3 3 10 20 30 40 50 60 70 80 90" | ./bin/blur_original
```
