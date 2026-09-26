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
| `tests/verify_vhdl.py` | Simula con GHDL e Icarus Verilog el VHDL/Verilog generado y los testbenches, y los compara bit a bit con el modelo en coma fija |
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

## Fase 5: Verilog y testbenches

* **Verilog-2001** del filtro FIR o IIR con los mismos puertos y la misma latencia que el
  VHDL (en pyfda solo existe la exportación Verilog de FIR con amaranth). Se comprueba
  con Icarus Verilog bit a bit frente al modelo, también con desbordamientos.
* **Testbenches autoverificables** en VHDL y Verilog (*Export HDL ... → with testbench*):
  impulso, escalones y ruido pseudoaleatorio, salida esperada calculada con el modelo en
  coma fija; imprimen `PASSED` o los errores muestra a muestra.
* El entero de cada coeficiente (formato `qint` de pyfda) ya aparece en la tabla de la
  pestaña Fixpoint junto a hex/bin/oct/CSD. pyfda solo tiene las estructuras FIR en
  forma directa e IIR DF1, las dos portadas.

## Fase 6: fórmula, fichero y espectrograma en y[n]

* Estímulo **Formula** como en pyfda (sintaxis de numexpr): índice `n`, tiempo `t = n / f_S`,
  parámetros A1, A2, f1, f2, phi1, phi2, T1, T2, N1, BW1, BW2, f_S, pi, e, operadores
  aritméticos, comparaciones y `& | ~`, y funciones (`sin`, `exp`, `where`, ...). Un
  evaluador propio en C++ sustituye a numexpr; los tests lo comparan con numexpr.
  Solo valores reales (pyfda también admite `j`).
* Estímulo **File**: carga CSV/texto, wav o npy (primera columna de datos que no sea el
  tiempo), normalizado opcionalmente a max |x| = 1 y escalado por A1. Con N = auto se usa
  la longitud del fichero.
* Pestaña **Spectrogram** (port de `scipy.signal.spectrogram`: PSD, magnitud o fase, en dB
  o lineal, NFFT y solapamiento, ventana de la pestaña Frequency), de x[n] o y[n], con
  barra de color. Comparado con scipy.

## Fase 7: ventanas, unidades y especificaciones

* Ventanas nuevas: Bartlett-Hann, Bohman, Cosine, Parzen, Triangular, **Dolph-Chebyshev**
  (atenuación de lóbulos laterales) y **DPSS/Slepian** (NW, autovector de la matriz
  tridiagonal por bisección e iteración inversa). Idénticas a `scipy.signal.get_window`.
  Sirven para el método de ventana, el análisis espectral y el espectrograma. Para el
  análisis espectral (como en pyfda) también **General Gaussian** con dos parámetros
  (forma p y σ en muestras).
* Unidades de frecuencia como en pyfda: f_S y **f_Ny** normalizadas, **mHz**, Hz, kHz, MHz y
  **GHz**. (La unidad `k` está desactivada en pyfda.)
* Especificaciones de amplitud en **dB, V o W** (`unit2lin`/`lin2unit` de pyfda, distintas
  para FIR e IIR en la banda de paso). Internamente y en el JSON se guardan en dB.

## Fase 8: información y visor de ventanas

* Pestaña **Info** (como `input_info` de pyfda): tipo, orden, nº de coeficientes y
  secciones, estabilidad (máx. |p|), fase mínima, |H(0)|, |H(f_S/2)|, máx. |H| y una
  tabla **especificaciones vs. conseguido** por banda (rizado en la banda de paso,
  atenuación en la de rechazo) marcando OK o fallo. Comprobado con `scipy.signal.freqz`.
* **Tools → Window viewer ...** (como `plot_fft_win` de pyfda): ventana en el tiempo,
  espectro con zero padding (lineal o dB) y sus propiedades: ganancia coherente, NENBW,
  pérdida de scalloping, ancho de banda a 3 y 6 dB y lóbulo lateral máximo. Comparado con
  los valores calculados en numpy.

## Fase 9: fichero de configuración

Como el `pyfda_user.conf` de pyfda, en un fichero INI legible
(`~/.config/pyfda/pyfda_cpp.ini`, en Windows `%APPDATA%\pyfda\pyfda_cpp.ini`), editable
también desde **File → Preferences ...**:

* **Sesión**: al cerrar se guardan el diseño (`pyfda_cpp_session.json`, junto al INI, con
  la coma fija), el tamaño de la ventana, los divisores y la pestaña; al arrancar se
  restauran (se puede desactivar).
* **Carpetas**: los diálogos recuerdan la última carpeta de filtros, de datos y de exportación.
* **Formato CSV** de las exportaciones (datos filtrados, análisis transitorio, tabla de
  coeficientes): separador `,`, `;` o tabulador y **coma decimal** (p. ej. para un Excel
  en español). Data Filt vuelve a leer esos ficheros.

`--config-dir carpeta` usa otra carpeta de configuración y `--quit` cierra la ventana al
momento (guardando la sesión), para las pruebas.

## Fase 10: señales complejas

Como en pyfda, el análisis transitorio y Data Filt trabajan con señales complejas:

* Estímulo **Exp (complex)**: A1·exp(j(2π f1 n + φ1)) + A2·exp(j(2π f2 n + φ2)).
* **A1, A2, DC y el ruido** admiten valores complejos (`1 - 3j`, `2j`, `exp(1j*pi/4)`); un
  ruido complejo añade ruido independiente a la parte real y a la imaginaria.
* **Fórmulas con `j`** (`exp(2j*pi*f1*n)`, `sqrt(-1+0j)`) y las funciones `real`, `imag`,
  `conj` y `complex`. Como numexpr, una fórmula real se sigue evaluando con números reales.
* En y[n] las partes **real e imaginaria** se dibujan en dos gráficas, el **espectro** y el
  **espectrograma** son de dos lados (-f_S/2 ... f_S/2) y el CSV exporta `x_re, x_im,
  y_re, y_im`. Los coeficientes son reales, así que la parte real y la imaginaria se
  filtran por separado (también en coma fija).
* **Datos complejos en ficheros**: arrays `.npy` complejos (complex64/128) y celdas CSV como
  `1+2j`, `(1+2j)` (formato de `numpy.savetxt`) o `1+2i` (Matlab). Data Filt dibuja la
  parte real (continua) y la imaginaria (discontinua), el espectro de dos lados y exporta
  las dos partes; el estímulo "File" también acepta datos complejos.

Comparado con numpy (fórmulas, estímulos), scipy (`spectrogram(..., return_onesided=False)`,
`sosfilt` de datos complejos) y la lectura de ficheros de numpy.

## Fase 11: gráfica 3D

Pestaña **3D** (como `plot_3d` de pyfda): |H(z)| sobre el plano z como **superficie**
coloreada (viridis, con sombreado) o **malla**, dentro del círculo unidad (rejilla polar) o
en |Re|, |Im| < 1,5, escala lineal o en dB con límites inferior y superior (la superficie
se recorta ahí), |H(f)| a lo largo del círculo unidad, el círculo unidad, polos (x) y ceros
(o) con sus tallos y barra de color. Se gira arrastrando con el ratón, se amplía con la
rueda y doble clic restablece la vista. |H(z)| está comparado con numpy.

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
`--save-filter f.json`, `--stim sine`, `--formula "sin(2*pi*f1*n)"`, `--stim-file x.csv`, `--tran-export t.csv`, `--fixpoint`,
`--export-hdl filtro.vhd|.v|.coe`, `--hdl-testbench`, `--data fichero`,
`--filter`, `--export salida.csv`, `--manual-ba "1,2,1/1,-0.5"`,
`--manual-zpk "0.5:0.5,0.5:-0.5/0.9/2"`, `--screenshot carpeta`, `--config-dir carpeta`, `--quit`.

## Pendiente para siguientes fases

* Edición de polos/ceros arrastrándolos en la gráfica P/Z.
* Traducciones.

## Licencias

pyfda: MIT. Las rutinas de diseño portadas de scipy.signal (incluido el algoritmo de
Remez de McClellan, Parks y Rabiner en la versión C de E. Kvaleberg) están bajo la
licencia BSD de SciPy.
