# pyfda con pestaña "Data Filt"

Copia de [pyfda](https://github.com/chipmuenk/pyfda) (commit `aa0eae30`) con una
pestaña nueva para aplicar el filtro diseñado a datos de un fichero.

## Uso

1. Arranca pyfda (`python -m pyfda.pyfdax` o `pyfdax.exe`).
2. Diseña el filtro en la parte izquierda (*Specs* → *DESIGN FILTER*).
3. Abre la pestaña **Data Filt** en la zona de gráficos.
4. **Load data ...**: elige un `.csv`, `.txt`, `.wav` o `.npy`.
   En `.csv` / `.txt` se detectan solos el separador (`;` `,` tabulador, espacios),
   la coma decimal, la cabecera y la codificación (UTF-8, UTF-8 de Excel, ANSI).
   Se saltan las líneas de comentario (`#`, `%`, `//`) y las de metadatos con otro
   número de campos, como las que añaden los osciloscopios.
5. **Column**: columna a filtrar. **Time**: columna de tiempo, o `n / f_S` para
   construir el eje con la frecuencia de muestreo del diseño.
6. **Filter data**: dibuja el dato filtrado encima del original.
   - *Zero phase*: filtra adelante y atrás (`filtfilt`), sin retardo.
   - *Spectrum*: añade los espectros de ambos.
   - *Export ...*: guarda tiempo, original y filtrado en un CSV.

Si cambias el diseño del filtro, el dato filtrado se recalcula al volver a la pestaña.
Las frecuencias del filtro se interpretan respecto a `f_S` del diseño, así que
ajusta `f_S` a la frecuencia de muestreo real de tus datos. Si eliges una columna
de tiempo (en segundos) y no coincide con `f_S`, aparece un aviso naranja con el
valor que hay que poner.

`datos_ejemplo.csv` contiene 10 Hz + 400 Hz + offset muestreados a 1 kHz.

## Ejecutable de Windows

El workflow `.github/workflows/build_pyfda_exe.yml` genera `pyfdax.exe` con
PyInstaller en cada push que toque `pyfda/` (o a mano desde la pestaña *Actions*).
Descárgalo desde el artefacto `pyfdax_win` de la ejecución.

Para generarlo en local en Windows:

```
cd pyfda
pip install -r requirements.txt pyinstaller pyinstaller-hooks-contrib
pip install . --no-deps
pyinstaller pyfdax.spec
```
