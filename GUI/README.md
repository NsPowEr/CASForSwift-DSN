# Detachable CAS GUI Lab

Questa cartella contiene il laboratorio GUI opzionale per test manuali del motore CAS.
Non e' parte del core matematico e non e' richiesta dal prodotto finale.

## Scopo

- Provare formule manualmente senza scrivere subito un test.
- Visualizzare output `Text`, `LaTeX`, `ASCII 2D` e campioni 2D.
- Restare completamente staccabile: rimuovere `GUI/` non deve rompere `src/`, `include/`, `test/` o `tools/cas_ui/`.

## Build

La GUI e' disabilitata di default:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build
```

Per abilitarla:

```bash
cmake -S . -B build-gui -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCAS_ENABLE_GUI=ON
cmake --build build-gui --target cas_gui
./build-gui/GUI/cas_gui
```

Serve Qt 6.7+ con `Quick` e `QuickControls2`. Se Qt non e' installato, il core resta compilabile perche questa cartella viene ignorata finche `CAS_ENABLE_GUI=OFF`.

## Confini

- La GUI linka `cas_core`; il core non deve mai includere header della GUI.
- Nessuna dipendenza Qt entra in `src/` o `include/cas/`.
- Nessuna Swift UI, FastAPI, DB, CloudKit, sync o AI client dentro questo repo.
- Plot e rendering avanzato devono passare da API/adattatori locali, non da hack nel symbolic core.

## Roadmap locale

1. Collegare input manuale a parse/simplify/format.
2. Aggiungere plot 2D client-side con `numeric::AdaptiveSampler`.
3. Aggiungere inspector AST read-only.
4. Aggiungere pannelli per `diff`, `integrate`, `limit`, `series`.
5. Solo dopo: valutare renderer LaTeX professionale e 3D.
