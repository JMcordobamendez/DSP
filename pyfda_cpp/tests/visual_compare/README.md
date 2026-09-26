# Comparación visual con pyfda

Abre el pyfda original y pyfda_cpp sin pantalla con los mismos seis diseños
(`designs.py`), guarda una captura de cada pestaña y genera parejas de imágenes
(pyfda a la izquierda, pyfda C++ a la derecha).

```
python3 -m venv venv && venv/bin/pip install numpy scipy matplotlib pyqt5 docutils mplcursors numexpr markdown pillow
venv/bin/pip install --no-deps -e <carpeta de pyfda>        # p. ej. la rama del PR #4
mkdir -p out/py
DATA_CSV=pyfda_cpp/examples/measurement.csv QT_QPA_PLATFORM=offscreen venv/bin/python shot_pyfda.py out/py
DATA_CSV=pyfda_cpp/examples/measurement.csv python3 shot_cpp.py <build>/pyfda_cpp out/cpp
venv/bin/python pair.py out                                 # -> out/pairs/<diseño>__<pestaña>.png
```

`shot_pyfda.py` guarda además `ba.json` y `shot_cpp.py` el filtro (`saved.json`) de
cada diseño, para comparar los coeficientes.
