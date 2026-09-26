# pyfda C++

Port a C++17 / Qt 6 de [pyfda](https://github.com/chipmuenk/pyfda) (Python Filter
Design Analysis Tool) incluyendo la pestaña **Data Filt** (filtrar una columna de un
CSV con el diseño actual). Es un port incremental: la aplicación ya permite diseñar
filtros IIR/FIR, ver sus respuestas, guardarlos y exportar los coeficientes, simular la
respuesta a distintos estímulos y filtrar datos medidos.

## Estructura

| Carpeta | Contenido |
|---|---|
| `src/core/` | Núcleo numérico sin dependencias (ni Qt): ports de `scipy.signal` y de la lógica de pyfda |
| `src/gui/` | Interfaz Qt 6 Widgets, gráficas propias con `QPainter` (sin matplotlib ni Qt Charts) |
| `tools/pyfda_cli.cpp` | CLI que expone el núcleo para compararlo con scipy |
| `tests/verify_scipy.py` | ~1400 comprobaciones contra scipy y contra el lector CSV en Python de pyfda |
| `examples/measurement.csv` | CSV de ejemplo (1 kHz, 10 Hz + 300 Hz + ruido) |
| `examples/lowpass_50Hz.json` | Filtro de ejemplo (elíptico, paso bajo 50 Hz a f_S = 1 kHz) |

## Qué está portado (fase 1)

* **Diseño IIR**: Butterworth, Chebyshev 1 y 2, elíptico y Bessel; paso bajo/alto/banda/
  banda eliminada; orden mínimo (`buttord`, `cheb1ord`, `cheb2ord`, `ellipord`, incluida
  la optimización `fminbound` para banda eliminada) u orden manual. Salida en zpk, b/a y
  secciones de segundo orden (`zpk2sos` con el emparejamiento `nearest` de scipy).
* **Diseño FIR**: método de ventana (`firwin`; rectangular, Bartlett, Hann, Hamming,
  Blackman, Blackman-Harris, Nuttall, flattop, Kaiser, gaussiana, Tukey; `kaiserord`) y
  equirizado Parks-McClellan (`remez`, portado del código C de scipy). Estimación del
  orden mínimo con `remezord` de pyfda (ichige, kaiser, herrmann).
* **Vistas**: |H(f)| (dB, V, W, con especificaciones sombreadas), fase, retardo de grupo,
  polos/ceros, respuesta al impulso/escalón, tabla de coeficientes (b/a, SOS, zpk; copiar
  y exportar CSV). Zoom con rueda o rectángulo, desplazamiento con Ctrl+arrastrar o botón
  central, doble clic para restablecer, menú contextual para copiar/guardar la imagen.
* **Data Filt** (port de `plot_data_filt.py`): carga csv/txt/wav/npy; el lector de texto
  detecta codificación (UTF-8, UTF-16, cp1252), delimitador, coma decimal, cabecera y
  líneas de comentario/metadatos igual que la versión Python; columna de tiempo o
  `n / f_S`; filtrado normal o de fase cero (`sosfilt`/`sosfiltfilt` para IIR,
  `lfilter`/`filtfilt` para FIR); espectro; aviso si la frecuencia de muestreo de los datos
  no coincide con `f_S`; exportación a CSV.

## Fase 2

* **Guardar y cargar filtros** (Archivo → Abrir / Guardar filtro, Ctrl+O / Ctrl+S) en un
  fichero JSON legible con las especificaciones, la unidad de frecuencia y los
  coeficientes. Al cargarlo se vuelve a diseñar desde las especificaciones y se avisa si
  los coeficientes guardados no coinciden. pyfda guarda en `.pkl`/`.npz` de Python, que no
  se pueden leer sin Python.
* **Exportar coeficientes** (Archivo → Exportar coeficientes, Ctrl+E, o el botón de la
  pestaña Coeffs) para MATLAB/Octave (`.m`), C (`.h`), Python/NumPy (`.py`) y CSV, con 17
  cifras significativas: los tests comprueban que cada formato se lee de vuelta bit a bit.
* **Análisis transitorio** (pestaña `y[n]`, port de `plot_tran.py` y
  `plot_tran_stim.py`): estímulos Dirac, sinc, gauss, rect, escalón (con error de
  asentamiento), seno, coseno, Dirichlet, chirp (lineal, cuadrático, logarítmico,
  hiperbólico), triángulo, diente de sierra, rectangular periódica, peine, AM, PM/FM y
  PWM, con versiones limitadas en banda (BL) como pyfda; ruido gaussiano, uniforme,
  enteros aleatorios, MLS y browniano; DC. Vista temporal (lineal o dB) y espectral con
  ventana (X, Y, potencia, |H| ideal, y respuesta en frecuencia a partir de un impulso).
  Exportación de n, t, x, y a CSV. Todos los estímulos deterministas coinciden con
  pyfda/scipy (error < 1e-12). Diferencias: la fase del chirp se interpreta en grados
  (pyfda pasa radianes a `scipy.signal.chirp`, que espera grados) y el ruido aleatorio no
  reproduce la secuencia de NumPy.

Diferencia con pyfda: en los FIR, `N` es siempre el **orden** (número de coeficientes − 1).
pyfda usa `N` como número de coeficientes en el método de ventana.

## Verificación

```
cmake -S . -B build -DPYFDA_BUILD_GUI=OFF && cmake --build build
python tests/verify_scipy.py build/pyfda_cli
```

Compara prototipos, polos/ceros, SOS, b/a, orden mínimo, ventanas, `firwin`, `remez`,
filtrado, respuesta en frecuencia, retardo de grupo, raíces, espectro, diseños completos
y la importación de CSV/wav/npy con scipy (tolerancia relativa ~1e-9). Si la variable
`PYFDA_DATA_FILT` apunta a `plot_data_filt.py`, también compara el lector CSV con la
función Python original.

## Compilar

Requiere CMake ≥ 3.16, un compilador C++17 y Qt 6 (solo el módulo Widgets).

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/pyfda_cpp
```

El workflow `.github/workflows/build_pyfda_cpp.yml` compila en Windows (MSVC + Qt 6.7),
pasa los tests contra scipy y sube `pyfda_cpp.exe` con las DLL de Qt como artefacto
`pyfda_cpp_win`.

Opciones de línea de comandos (para pruebas): `--load-filter f.json`,
`--save-filter f.json`, `--stim sine`, `--tran-export t.csv`, `--data fichero`,
`--filter`, `--export salida.csv`, `--screenshot carpeta`.

## Pendiente para siguientes fases

* En el análisis transitorio: estímulos complejos (exp), fórmula libre, espectrograma,
  datos de fichero como estímulo y condiciones iniciales. Gráfica 3D.
* Coma fija: cuantización, `fixpoint_widgets`, simulación y exportación HDL (amaranth).
* Edición manual de coeficientes y de polos/ceros, métodos *Moving Average*, *Delay* y
  *Manual*.
* Exportar coeficientes cuantizados (COE de Xilinx, VHDL, ...), junto con la coma fija.
* Visor de ventanas FFT y el resto de ventanas de pyfda, más unidades (f_Ny, k),
  especificaciones de amplitud en V/W, pestaña de información, fichero de configuración
  y traducciones.
* Datos complejos en Data Filt (ahora solo reales).

## Licencias

pyfda: MIT. Las rutinas de diseño portadas de scipy.signal (incluido el algoritmo de
Remez de McClellan, Parks y Rabiner en la versión C de E. Kvaleberg) están bajo la
licencia BSD de SciPy.
