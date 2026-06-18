# 🔍 OBIETTIVO: Audit di Parità, Integrità e Potenza Matematica (Fasi F0-F7.5)

Sei un **Lead CAS Architect e Principal Auditor**. Il tuo compito è **certificare la reale potenza matematica** del `REAL CAS ENGINE C++` rispetto al target di parità (HP Prime G2 / Giac).

**È vietata l`accettazione passiva dei ledger o dei report storici.** Il tuo audit deve essere **matematicamente empirico, spietato e profondo**. Qualsiasi costrutto che appaia come "finito" ma celi un hardcode, una limitazione di dominio non documentata, un fallimento silenzioso (error masking) o un bug architetturale, deve essere denunciato.

---

## 🛠️ DIRETTIVE DI ISPEZIONE (Il Metodo)

1. **Uso pervasivo del Grafo della Conoscenza (`graphify-out/`)**:
   - Inizia sempre eseguendo query sul grafo per mappare le dipendenze occulte:
     - `graphify query "symbolic simplification cycle"`
     - `graphify explain "assumption propagation"`
   - Usa il grafo per stanare codice morto, cicli di dipendenza o funzioni "god object".

2. **Orchestrazione Avanzata di Subagenti (Deep Dive)**:
   - **Non lavorare da solo.** Devi orchestrare indagini parallele utilizzando la tool `invoke_agent`.
   - Delega esplorazioni mirate ai subagenti specializzati (es. `codebase_investigator`, `cas-calculus-dev`, `cas-algebra-dev`, `cas-symbolic-dev`).
   - Sfrutta i subagenti per ispezionare gli angoli più nascosti del codice, analizzare l`implementazione degli algoritmi e individuare debiti che sfuggirebbero a un`analisi di superficie.
   - Consolida tutti i responsi in un **Singolo Report di Stato unificato, visualmente appetibile e dettagliato**.

3. **Analisi Spietata del Gap (La via verso il 100%)**:
   - Per **ogni singola fase o modulo (F1-F7)** che non ha raggiunto il 100% o presenta debiti, devi produrre una diagnosi profonda:
     - **La Motivazione:** Perché la fase non è al 100%? Quale ostacolo matematico, prestazionale o architetturale lo impedisce?
     - **Le Dipendenze Mancanti:** Quali strutture dati, astrazioni o algoritmi collaterali servono per sbloccare la fase?
     - **Il Piano d`Azione:** Qual è l`esatta sequenza di azioni (algoritmi, refactor, implementazioni) necessaria per chiudere il gap e arrivare al 100% di potenza reale?

4. **Caccia al Falso Positivo e Robustezza Strutturale**:
   - Ispeziona le euristiche (es. GCD euristico) e i _bail-out_ (`CASErrorKind::Undefined`, `nullopt`).
   - Verifica i test: coprono l`intero dominio (C, estensioni algebriche, singolarità) o solo polinomi in Q[x]?
   - Controlla lo stato reale del waiver sui file monolitici (>500 LOC).

---

## 📊 STATO ATTUALE: MAPPA DI PARITÀ (Pre-Audit)

> _Nota per l`Auditor: Questa è la dichiarazione attuale del progetto. Demoliscila, verificala e sostituiscila con il tuo Master Report (indicando in % lo stato reale dove 100% significa totalmente privo di debiti)._

| Modulo / Dominio | Completamento Stimato | Potenza Reale vs HP Prime | Debiti Aperti (Da verificare) |
| :--- | :---: | :--- | :--- |
| **F1: Foundation** | 🟢 95% | **Alta**. Toom-3, Karatsuba, Immutable AST. | HPP-F1.1-MUL (Schönhage-Strassen), `is_neg_infinity` legacy form. |
| **F2: Univariate Alg.**| 🟢 98% | **Molto Alta**. Berlekamp, Hensel. | Nessuno noto. |
| **F3: Multivar & Ext.**| 🟡 80% | **Media**. Brown GCD, Trager primitive, FGLM. | Wang factorization, Galois deg >= 6, iterazioni torri complesse. |
| **F4: Linear Algebra** | 🟢 90% | **Alta**. Bareiss, MGS QR, Cholesky, Smith. | HPP-F4.1-QR-HOUSEHOLDER (Householder QR simbolico). |
| **F5: Calculus**       | 🟡 75% | **Media**. Derivate, limiti Gruntz. | Risch Hermite reduction completa, Risch transcendental, Equazioni differenziali. |
| **F6: Numeric, CAD**   | 🟠 50% | **Bassa**. MPFR, branch cuts parziali. | CAD generale, Multi-sheet Riemann. |
| **F7: Acceptance**     | 🟢 94.5%| **Aggregato Alto**. 94.5% non-skip. | Hypergeometric `_pF_q` completo. |

**Architettura & Code Quality:**
- **Anti-Monolith**: 🔴 Fallito. Esiste `HC-F8-MONOLITH-WAIVER` (27 file >500 LOC).

---

## 🚀 I PROSSIMI PASSI E L`ORIZZONTE FUTURO

Il tuo output finale deve includere un **Master Action Plan** diviso in tre sezioni:

### PARTE A: Remediazione Immediata e Gap Analysis (Risultati dell`Audit)
Sintetizza i risultati dei subagenti in un singolo report a colpo d`occhio. Correggi le percentuali. Spiega **esattamente cosa manca, perché manca e come implementarlo** per portare ogni fase storica al 100%.

### PARTE B: Roadmap Fase 8 (Il Cammino verso la Parità HP Prime)
Prioritizza le lavorazioni strategiche a breve/medio termine:
1. **F8.0 Prerequisito Tecnico**: Esecuzione del piano `ANTI_MONOLITHIC_REPORT.md` per i 14 file Tier-1 (>600 LOC).
2. **F8.1 Integrazione Simbolica**: Chiusura del Bronstein gap (Hermite reduction, trascendentale completo).
3. **F8.2 Algebra Avanzata**: Implementazione CAD (Cylindrical Algebraic Decomposition), Galois >= 6, riconoscimento `_pF_q`.
4. **F8.3 Analisi Complessa**: Gestione branch cuts rigorosa e superfici di Riemann multi-foglio.

### PARTE C: L`Orizzonte Post-Fase 8 (La Visione Definitiva per l`Equivalenza Totale)
**Ragiona profondamente:** Una volta conclusa la Fase 8, cosa mancherà per poter affermare in modo inoppugnabile che il nostro CAS C++ eguaglia o supera in toto il motore Giac dell`HP Prime G2?
- Individua le sfide ingegneristiche finali (es. JIT Compilation per le espressioni valutate di frequente, ottimizzazioni assembly low-level per le code di BigInt, sistemi di cache distribuiti avanzati o pattern matching guidato da euristiche probabilistiche).
- Delinea le eventuali mancanze che potrebbero sorgere in futuro nel codice attualmente sviluppato, delineando il "Piano Finale" (Fase 9) per certificare l`equivalenza totale a livello di sistema, usabilità e performance.

---
**Comando per iniziare:** Attiva i tuoi subagenti, esplora il grafo e il codice in ogni angolo, definisci i percorsi per il 100% ed elabora il Report Finale di Stato.
