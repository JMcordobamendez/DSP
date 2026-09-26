# pyfda C++

Port a C++17 / Qt 6 de [pyfda](https://github.com/chipmuenk/pyfda) (Python Filter
Design Analysis Tool) incluyendo la pestaña **Data Filt** (filtrar una columna de un
CSV con el diseño actual). Esta es la **fase 1** de un port incremental: la aplicación
ya funciona para diseñar filtros IIR/FIR, ver sus respuestas y filtrar datos medidos.

## Estructura

| Carpeta | Contenido |
|---|---|
| `src/core/` | Núcleo numérico sin dependencias (ni Qt): ports de `scipy.signal` y de la lógica de pyfda |
| `src/gui/` | Interfaz Qt 6 Widgets, gráficas propias con `QPainter` (sin matplotlib ni Qt Charts) |
| `tools/pyfda_cli.cpp` | CLI que expone el núcleo para compararlo con scipy |
| `tests/verify_scipy.py` | ~1300 comprobaciones contra scipy y contra el lector CSV en Python de pyfda |
| `examples/measurement.csv` | CSV de ejemplo (1 kHz, 10 Hz + 300 Hz + ruido) |

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

Opciones de línea de comandos (para pruebas): `--data fichero`, `--filter`,
`--export salida.csv`, `--screenshot carpeta`.

## Pendiente para siguientes fases

* Pestaña de análisis transitorio (`plot_tran`: estímulos, ventanas, FFT) y gráfica 3D.
* Coma fija: cuantización, `fixpoint_widgets`, simulación y exportación HDL (amaranth).
* Edición manual de coeficientes y de polos/ceros, métodos *Moving Average*, *Delay* y
  *Manual*.
* Guardar/cargar diseños y exportar coeficientes a otros formatos (COE, VHDL, ...).
* Visor de ventanas FFT y el resto de ventanas de pyfda, más unidades (f_Ny, k),
  especificaciones de amplitud en V/W, pestaña de información, fichero de configuración
  y traducciones.
* Datos complejos en Data Filt (ahora solo reales).

## Licencias

pyfda: MIT. Las rutinas de diseño portadas de scipy.signal (incluido el algoritmo de
Remez de McClellan, Parks y Rabiner en la versión C de E. Kvaleberg) están bajo la
licencia BSD de SciPy.
