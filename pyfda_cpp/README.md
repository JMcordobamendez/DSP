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
| `tests/verify_scipy.py` | ~1700 comprobaciones contra scipy/pyfda y contra el lector CSV en Python de pyfda |
| `tests/verify_vhdl.py` | Simula con GHDL el VHDL generado y lo compara bit a bit con el modelo en coma fija |
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

## Fase 3: coma fija

* **Cuantizador** (port de `Fixed` de `pyfda_fix_lib.py`): formatos Q con WI bits enteros,
  WF fraccionarios y signo; cuantización floor, round (a par, como `numpy.round`), fix,
  ceil o ninguna; desbordamiento por saturación, *wrap-around* en complemento a dos o
  ninguno, con contadores. Coincide bit a bit con pyfda (comparado con 285 formatos).
* **Pestaña Fixpoint**: formatos de entrada, coeficientes b y a, acumulador y salida, con
  modo automático (bits enteros de los coeficientes y acumulador sin redondeo de los
  productos con bits de guarda); tabla de coeficientes cuantizados (valor, error, entero
  y hex/bin/oct/CSD) y |H(f)| ideal frente a la de los coeficientes cuantizados.
* **Simulación en coma fija** en la pestaña `y[n]`: respuesta en coma fija, respuesta en
  coma flotante como referencia, desbordamientos por etapa y error máximo. FIR en forma
  directa (igual que `fir_df_pyfixp`, bit a bit); IIR como cascada de secciones de segundo
  orden en forma directa 1 (pyfda usa `iir_df1` con la función de transferencia completa,
  que se vuelve inestable en órdenes altos).
* **Exportación**: fichero COE de Xilinx (FIR Compiler) y **VHDL sintetizable**
  (`numeric_std`, VHDL-2008) del filtro FIR o IIR. `tests/verify_vhdl.py` lo simula con
  GHDL y comprueba que la salida es idéntica bit a bit al modelo, también con
  desbordamientos. pyfda genera Verilog con amaranth (Python); aquí se genera VHDL sin
  dependencias.
* Los ajustes de coma fija se guardan en el fichero JSON del filtro.

Diferencia con pyfda: en los FIR, `N` es siempre el **orden** (número de coeficientes − 1).
pyfda usa `N` como número de coeficientes en el método de ventana.

## Fase 4: edición manual, Moving Average, Delay y Manual

* **Moving Average** (port de `ma.py`): paso bajo, paso alto, paso banda y elimina banda
  con varias etapas en cascada y normalización opcional; orden mínimo para paso bajo y
  paso alto a partir de F_SB y A_SB. Los coeficientes y los ceros coinciden con el código
  de pyfda (`calc_ma`). Dos diferencias: en paso banda y elimina banda los ceros se
  calculan como raíces de b (los de pyfda no corresponden a sus coeficientes) y
  `ceil_odd` devuelve el impar ≥ x como dice su documentación (el de pyfda suma 2 a los
  impares y el orden crece en cada rediseño).
* **Delay**: N retardos, H(z) = z<sup>-N</sup>.
* **Manual** (FIR o IIR): en la pestaña **Coeffs**, *Edit* permite escribir b y a o los
  ceros, polos y la ganancia k, añadir y borrar filas y *Apply* (Ctrl+Enter) diseña un
  filtro *Manual* con ellos. El tipo FIR/IIR se decide por los coeficientes, a se
  normaliza a a[0] = 1 y si hay más ceros que polos se añaden polos en el origen. Los
  polos/ceros complejos deben ir en pares conjugados. Al elegir *Manual* en el menú se
  parte de los coeficientes del filtro actual, como en pyfda.
* Los filtros manuales se guardan en el JSON con sus coeficientes o, si se introdujeron
  como polos/ceros, con `"zpk"` para no perder precisión.

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
`--save-filter f.json`, `--stim sine`, `--tran-export t.csv`, `--fixpoint`,
`--export-hdl filtro.vhd|.coe`, `--data fichero`,
`--filter`, `--export salida.csv`, `--manual-ba "1,2,1/1,-0.5"`,
`--manual-zpk "0.5:0.5,0.5:-0.5/0.9/2"`, `--screenshot carpeta`.

## Pendiente para siguientes fases

* En el análisis transitorio: estímulos complejos (exp), fórmula libre, espectrograma,
  datos de fichero como estímulo y condiciones iniciales. Gráfica 3D.
* Coma fija: formato de visualización entero (`qint`), otras estructuras (FIR transpuesta,
  IIR DF2), exportación Verilog y testbench VHDL.
* Edición de polos/ceros arrastrándolos en la gráfica P/Z.
* Visor de ventanas FFT y el resto de ventanas de pyfda, más unidades (f_Ny, k),
  especificaciones de amplitud en V/W, pestaña de información, fichero de configuración
  y traducciones.
* Datos complejos en Data Filt (ahora solo reales).

## Licencias

pyfda: MIT. Las rutinas de diseño portadas de scipy.signal (incluido el algoritmo de
Remez de McClellan, Parks y Rabiner en la versión C de E. Kvaleberg) están bajo la
licencia BSD de SciPy.
