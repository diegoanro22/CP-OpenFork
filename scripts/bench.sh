#!/usr/bin/env bash
#
# Mide tiempo, speedup y eficiencia. Cada quien lo corre en su maquina y el
# resultado queda en docs/resultados/<algoritmo>_<nombre>.csv.
#
# Uso: ./scripts/bench.sh [nombre]
# Para acortarlo: REPS=1 HILOS="1 2 4" N_MAT=512 M_IMG=1920 N_IMG=1080 ./scripts/bench.sh

set -euo pipefail
cd "$(dirname "$0")/.."

USUARIO="${1:-$(whoami)}"
REPS="${REPS:-3}"
N_MAT="${N_MAT:-1024}"
M_IMG="${M_IMG:-7680}"
N_IMG="${N_IMG:-4320}"

# Potencias de 2 hasta los hilos que tenga la maquina.
NUCLEOS=$(getconf _NPROCESSORS_ONLN)
if [ -z "${HILOS:-}" ]; then
    HILOS=""
    p=1
    while [ "$p" -le "$NUCLEOS" ]; do
        HILOS="$HILOS $p"
        p=$((p * 2))
    done
fi

make --no-print-directory all

mkdir -p docs/resultados

# Datos de la maquina, para el informe.
INFO="docs/resultados/maquina_${USUARIO}.txt"
{
    echo "Integrante : $USUARIO"
    echo "Fecha      : $(date '+%Y-%m-%d %H:%M')"
    echo "Nucleos    : $NUCLEOS hilos logicos"
    echo "CPU        : $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ //' || uname -m)"
    echo "Compilador : $(gcc --version | head -1)"
    echo "SO         : $(uname -sr)"
} > "$INFO"

mediana() {
    sort -g | awk '{v[NR]=$1} END{ print (NR%2) ? v[(NR+1)/2] : (v[NR/2]+v[NR/2+1])/2 }'
}

# medir <hilos> <comando...> -> "<mediana> <checksum>"
medir() {
    local hilos="$1"; shift
    local tiempos="" linea chk=""
    for _ in $(seq 1 "$REPS"); do
        linea=$(OMP_NUM_THREADS="$hilos" "$@")
        tiempos="$tiempos$(sed -n 's/.*tiempo=\([0-9.]*\).*/\1/p' <<<"$linea")"$'\n'
        chk=$(sed -n 's/.*checksum=\([0-9]*\).*/\1/p' <<<"$linea")
    done
    echo "$(printf '%s' "$tiempos" | mediana) $chk"
}

# Dos speedups: total contra el secuencial original, y paralelo contra la misma
# variante con 1 hilo (ese aisla lo que aporto OpenMP).
barrido() {
    local etiqueta="$1" tamano="$2" csv="$3" binario="$4"; shift 4
    local variantes=("$@")

    echo
    echo "=== $etiqueta ($tamano), $REPS repeticiones, mediana ==="
    printf '%-9s %-5s %11s  %9s %7s  %9s %7s\n' \
        variante hilos tiempo sp_total ef_total sp_par ef_par
    echo "algoritmo,variante,tamano,hilos,tiempo_s,speedup_total,eficiencia_total,speedup_paralelo,eficiencia_paralela" > "$csv"

    read -r t_base chk_base < <(medir 1 "$binario" seq "${ARGS[@]}")
    printf '%-9s %-5s %9.4f s  %8.2fx %6.1f%%  %8s %7s\n' seq 1 "$t_base" 1 100 "-" "-"
    echo "$etiqueta,seq,$tamano,1,$t_base,1.0000,1.0000,,"  >> "$csv"

    for variante in "${variantes[@]}"; do
        local t_un_hilo=""
        for h in $HILOS; do
            read -r t chk < <(medir "$h" "$binario" "$variante" "${ARGS[@]}")

            if [ "$chk" != "$chk_base" ]; then
                echo "ERROR: checksum distinto en $variante con $h hilos" >&2
                echo "       esperado $chk_base, obtenido $chk" >&2
                exit 1
            fi

            # La primera pasada siempre es con 1 hilo.
            [ -z "$t_un_hilo" ] && t_un_hilo="$t"

            read -r spt eft spp efp < <(awk -v b="$t_base" -v u="$t_un_hilo" -v t="$t" -v p="$h" \
                'BEGIN{ st=b/t; sp=u/t; printf "%.4f %.4f %.4f %.4f\n", st, st/p, sp, sp/p }')

            printf '%-9s %-5s %9.4f s  %8.2fx %6.1f%%  %8.2fx %6.1f%%\n' \
                "$variante" "$h" "$t" "$spt" \
                "$(awk -v e="$eft" 'BEGIN{print e*100}')" "$spp" \
                "$(awk -v e="$efp" 'BEGIN{print e*100}')"
            echo "$etiqueta,$variante,$tamano,$h,$t,$spt,$eft,$spp,$efp" >> "$csv"
        done
    done

    echo "-> $csv"
}

ARGS=("$N_MAT")
barrido matrices "N=$N_MAT" "docs/resultados/matrices_${USUARIO}.csv" \
    ./bin/matrices ijk ikj tiled

ARGS=("$M_IMG" "$N_IMG")
barrido blur "${M_IMG}x${N_IMG}" "docs/resultados/blur_${USUARIO}.csv" \
    ./bin/blur filas collapse dynamic

echo
echo "Listo. Checksums verificados contra la version secuencial."
echo "Datos de la maquina en $INFO"
echo "Para las graficas: python3 scripts/graficar.py $USUARIO"
