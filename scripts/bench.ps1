# Mide tiempo, speedup y eficiencia en Windows.
# Uso: .\scripts\bench.ps1 [nombre]
# O con variables: $env:REPS="1"; $env:HILOS="1 2 4"; .\scripts\bench.ps1

$ErrorActionPreference = "Stop"

# Cambiar al directorio raíz del proyecto
Set-Location (Join-Path $PSScriptRoot "..")

$USUARIO = if ($args[0]) { $args[0] } else { $env:USERNAME }
$REPS = if ($env:REPS) { [int]$env:REPS } else { 3 }
$N_MAT = if ($env:N_MAT) { $env:N_MAT } else { "1024" }
$M_IMG = if ($env:M_IMG) { $env:M_IMG } else { "7680" }
$N_IMG = if ($env:N_IMG) { $env:N_IMG } else { "4320" }

# Obtener número de hilos lógicos
$NUCLEOS = (Get-CimInstance Win32_ComputerSystem).NumberOfLogicalProcessors

if (-not $env:HILOS) {
    $HILOS_LIST = @()
    $p = 1
    while ($p -le $NUCLEOS) {
        $HILOS_LIST += $p
        $p = $p * 2
    }
} else {
    $HILOS_LIST = $env:HILOS -split '\s+' | ForEach-Object { [int]$_ }
}

# Compilar proyecto usando mingw32-make o make
if (Get-Command mingw32-make -ErrorAction SilentlyContinue) {
    mingw32-make --no-print-directory all
} else {
    make --no-print-directory all
}

New-Item -ItemType Directory -Force -Path "docs/resultados" | Out-Null

# Información de la máquina
$INFO = "docs/resultados/maquina_${USUARIO}.txt"
$cpuInfo = (Get-CimInstance Win32_Processor | Select-Object -First 1).Name
$gccVersion = (gcc --version | Select-Object -First 1)
$osVersion = (Get-CimInstance Win32_OperatingSystem).Caption

@"
Integrante : $USUARIO
Fecha      : $(Get-Date -Format 'yyyy-MM-dd HH:mm')
Nucleos    : $NUCLEOS hilos logicos
CPU        : $cpuInfo
Compilador : $gccVersion
SO         : $osVersion
"@ | Out-File -FilePath $INFO -Encoding utf8

function Obtener-Mediana ($lista) {
    $ordenados = $lista | Sort-Object
    $count = $ordenados.Count
    if ($count % 2 -ne 0) {
        return $ordenados[[math]::Floor($count / 2)]
    } else {
        return ($ordenados[$count / 2 - 1] + $ordenados[$count / 2]) / 2
    }
}

function Medir ($hilos, $exe, $variante, $argumentos) {
    $tiempos = @()
    $chk = ""
    
    $env:OMP_NUM_THREADS = $hilos
    
    for ($i = 0; $i -lt $REPS; $i++) {
        $salida = & $exe $variante $argumentos 2>&1 | Out-String
        
        if ($salida -match 'tiempo=([0-9.]+)') {
            $tiempos += [double]$Matches[1]
        }
        if ($salida -match 'checksum=([0-9]+)') {
            $chk = $Matches[1]
        }
    }
    
    $med = Obtener-Mediana $tiempos
    return @($med, $chk)
}

function Barrido ($etiqueta, $tamano, $csv, $binario, $variantes, $cmdArgs) {
    Write-Host "`n=== $etiqueta ($tamano), $REPS repeticiones, mediana ==="
    "{0,-9} {1,-5} {2,11}  {3,9} {4,7}  {5,9} {6,7}" -f "variante", "hilos", "tiempo", "sp_total", "ef_total", "sp_par", "ef_par" | Write-Host
    "algoritmo,variante,tamano,hilos,tiempo_s,speedup_total,eficiencia_total,speedup_paralelo,eficiencia_paralela" | Out-File -FilePath $csv -Encoding utf8

    $resBase = Medir 1 $binario "seq" $cmdArgs
    $t_base = [double]$resBase[0]
    $chk_base = $resBase[1]

    "{0,-9} {1,-5} {2,9:F4} s  {3,8:F2}x {4,6:F1}%  {5,8} {6,7}" -f "seq", 1, $t_base, 1.0, 100.0, "-", "-" | Write-Host
    "$etiqueta,seq,$tamano,1,$t_base,1.0000,1.0000,," | Out-File -FilePath $csv -Append -Encoding utf8

    foreach ($v in $variantes) {
        $t_un_hilo = $null
        foreach ($h in $HILOS_LIST) {
            $res = Medir $h $binario $v $cmdArgs
            $t = [double]$res[0]
            $chk = $res[1]

            if ($chk -ne $chk_base) {
                Write-Error "ERROR: checksum distinto en $v con $h hilos`n       esperado $chk_base, obtenido $chk"
                exit 1
            }

            if ($null -eq $t_un_hilo) { $t_un_hilo = $t }

            $spt = $t_base / $t
            $eft = $spt / $h
            $spp = $t_un_hilo / $t
            $efp = $spp / $h

            "{0,-9} {1,-5} {2,9:F4} s  {3,8:F2}x {4,6:F1}%  {5,8:F2}x {6,6:F1}%" -f $v, $h, $t, $spt, ($eft * 100), $spp, ($efp * 100) | Write-Host
            "$etiqueta,$v,$tamano,$h,$t,$spt,$eft,$spp,$efp" | Out-File -FilePath $csv -Append -Encoding utf8
        }
    }

    Write-Host "-> $csv"
}

# Ejecutar pruebas
$binMatrices = ".\bin\matrices.exe"
Barrido "matrices" "N=$N_MAT" "docs/resultados/matrices_${USUARIO}.csv" $binMatrices @("ijk", "ikj", "tiled") @($N_MAT)

$binBlur = ".\bin\blur.exe"
Barrido "blur" "${M_IMG}x${N_IMG}" "docs/resultados/blur_${USUARIO}.csv" $binBlur @("filas", "collapse", "dynamic") @($M_IMG, $N_IMG)

Write-Host "`nListo. Checksums verificados contra la version secuencial."
Write-Host "Datos de la maquina en $INFO"
Write-Host "Para las graficas: python scripts/graficar.py $USUARIO"