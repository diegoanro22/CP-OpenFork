#!/usr/bin/env python3
"""Graficas de speedup y eficiencia a partir de los CSV del benchmark.

Uso: python3 scripts/graficar.py [nombre]
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = Path(__file__).resolve().parent.parent
USUARIO = sys.argv[1] if len(sys.argv) > 1 else "diego"


def leer(csv_path):
    """{variante: [(hilos, tiempo, sp_total, sp_par, ef_par), ...]}"""
    series = defaultdict(list)
    with open(csv_path) as f:
        for fila in csv.DictReader(f):
            if fila["variante"] == "seq":
                continue
            series[fila["variante"]].append((
                int(fila["hilos"]),
                float(fila["tiempo_s"]),
                float(fila["speedup_total"]),
                float(fila["speedup_paralelo"]),
                float(fila["eficiencia_paralela"]) * 100,
            ))
    for v in series.values():
        v.sort()
    return series


def graficar(csv_path, titulo, salida):
    series = leer(csv_path)
    if not series:
        return False

    fig, (izq, der) = plt.subplots(1, 2, figsize=(12, 4.5))
    hilos = [p[0] for p in next(iter(series.values()))]

    izq.plot(hilos, hilos, "--", color="gray", linewidth=1, label="ideal (lineal)")
    for variante, puntos in series.items():
        izq.plot([p[0] for p in puntos], [p[3] for p in puntos], "o-", label=variante)
    izq.set_xscale("log", base=2)
    izq.set_yscale("log", base=2)
    izq.set_xticks(hilos)
    izq.set_xticklabels(hilos)
    izq.set_xlabel("Hilos")
    izq.set_ylabel("Speedup vs. la misma variante con 1 hilo")
    izq.set_title("Speedup aportado por OpenMP")
    izq.grid(True, alpha=0.3)
    izq.legend()

    der.axhline(100, ls="--", color="gray", linewidth=1)
    for variante, puntos in series.items():
        der.plot([p[0] for p in puntos], [p[4] for p in puntos], "o-", label=variante)
    der.set_xscale("log", base=2)
    der.set_xticks(hilos)
    der.set_xticklabels(hilos)
    der.set_ylim(0, 140)
    der.set_xlabel("Hilos")
    der.set_ylabel("Eficiencia (%)")
    der.set_title("Eficiencia = speedup / hilos")
    der.grid(True, alpha=0.3)
    der.legend()

    fig.suptitle(titulo)
    fig.tight_layout()
    fig.savefig(salida, dpi=130)
    plt.close(fig)
    print(f"-> {salida}")
    return True


res = RAIZ / "docs" / "resultados"
gra = RAIZ / "docs" / "graficas"
gra.mkdir(parents=True, exist_ok=True)

hecho = False
for algo, titulo in [("matrices", "Multiplicacion de matrices (N=1024)"),
                     ("blur", "Filtro de desenfoque (7680x4320)")]:
    origen = res / f"{algo}_{USUARIO}.csv"
    if origen.exists():
        hecho |= graficar(origen, f"{titulo} - {USUARIO}", gra / f"{algo}_{USUARIO}.png")
    else:
        print(f"(no existe {origen}, se omite)")

if not hecho:
    sys.exit("No se genero ninguna grafica. Corra antes ./scripts/bench.sh")
