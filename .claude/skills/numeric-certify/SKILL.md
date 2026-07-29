---
name: numeric-certify
description: Certificato numerico multi-punto ad alta precisione (mpmath) per validare un risultato matematico del CAS — identità, integrale (D(F)=f), limite, serie, fattorizzazione. Metodo provato nel progetto (spec A7 verificata mpmath; certificati multi-punto in memoria). Non sostituisce i test strutturali: li integra quando l'equivalenza simbolica è dubbia.
---

# numeric-certify

Valida un claim matematico con campionamento numerico ad alta precisione
(mpmath, precisione e punti configurabili). Complementare al confronto
strutturale: da usare quando `simplify(a-b)==0` non chiude (memoria:
`polynomial_normal_form` non idempotente sotto `op_active`; simplify non fonde
`e^a·e^b`) o per validare formule di una spec prima di implementarle.

## Quando usare

- Verificare un fix/feature: risultato CAS vs forma attesa.
- Validare formule in una bozza di spec (obbligatorio per `spec-researcher`).
- Arbitrare `ORACLE_DIVERGE` (Maxima ≠ Giac ≠ CAS): il numero decide.
- Smascherare silent-wrong: identità che "sembrano" giuste.

## Protocollo

1. **Formula del certificato**: riduci il claim a `LHS(x…) − RHS(x…) ≈ 0`.
   - Integrale indefinito: certifica `d/dx F(x) − f(x)` (derivata numerica o
     simbolica mpmath — MAI fidarsi solo di F "che assomiglia").
   - Integrale definito: `mpmath.quad` vs valore chiuso.
   - Limite: campiona avvicinamento (sequenze geometriche verso il punto).
   - Serie/somma: `mpmath.nsum` vs forma chiusa.

2. **Punti di campionamento** (anti falso-positivo):
   - ≥ 7 punti, IRRAZIONALI e trascendenti misti (es. `sqrt(2)+1/3`, `pi/7`,
     `e/3`, `1+sqrt(5)`), MAI solo interi piccoli (0,1,2 nascondono
     cancellazioni accidentali).
   - Includi punti vicini a singolarità/branch cut dichiarati, da entrambi i
     lati se il dominio lo consente.
   - Multi-variabile: griglia casuale seedata dall'hash dell'input (Cat-6:
     mai seed fisso arbitrario).

3. **Precisione**: `mp.dps = 50` default; raddoppia (100) e ripeti — un vero
   zero scala, un falso zero da cancellazione no. Soglia: `|err| < 10^(-dps/2)`.

4. **Branch cuts**: per log/radici/potenze complesse, dichiara la convenzione
   (principal branch) e campiona su C, non solo su R, se il claim è complesso.

5. **Verdetto**:
   - `CERTIFIED`: tutti i punti sotto soglia a entrambe le precisioni.
   - `REFUTED`: un punto sopra soglia stabile → il claim è FALSO, riporta il
     controesempio esatto (punto + valori LHS/RHS).
   - `INCONCLUSIVE`: instabilità numerica/singolarità — restringi dominio o
     alza precisione; MAI spacciare inconclusive per certified.

## Esempio scheletro

```python
from mpmath import mp, mpf, sqrt, pi, e, log, exp, quad, diff
mp.dps = 50
pts = [sqrt(2)+mpf(1)/3, pi/7, e/3, 1+sqrt(5), mpf(7)/11, pi*e/5, sqrt(3)/2]
def LHS(x): ...
def RHS(x): ...
worst = max(abs(LHS(x)-RHS(x)) for x in pts)
assert worst < mpf(10)**(-mp.dps//2), f"REFUTED: worst={worst}"
```

## Regole

- Un certificato numerico NON è una dimostrazione: `CERTIFIED` = "nessun
  controesempio trovato". Per chiudere task serve comunque il criterio della
  spec (test strutturali/gate).
- `REFUTED` batte qualunque test verde: apri/riapri la task, non "aggiustare"
  i punti di campionamento (REGOLA 0.2).
- Script temporanei in `$CLAUDE_JOB_DIR/tmp` o `/tmp`, mai committati.
