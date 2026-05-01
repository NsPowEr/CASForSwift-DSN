# CAS Engine — Costituzione Tecnica per Agenti e Sviluppatori

> Questo file è la **Legge Suprema** del progetto. Ogni agente DEVE leggerlo e rispettarlo.
> Ignorare queste regole porterà al rigetto immediato del codice (Architectural Fail).

---

## REGOLE ARCHITETTURALI NON NEGOZIABILI

### 1. Aritmetica Esatta e BigInt a Limbs
- È **vietato** l'uso di `int64_t` o `double` per il calcolo simbolico.
- Il tipo `BigInt` deve usare un'implementazione a **Limbs** (`std::vector<uint32_t>` o `uint64_t`) con aritmetica bit-a-bit. 
- È **vietata** l'aritmetica basata su stringhe decimali (troppo lenta, $O(N^2)$).

### 2. Structural Sharing (Regola d'Oro delle Performance)
- L'AST è immutabile. Se una funzione (es: `simplify`, `substitute`) non modifica un nodo, **DEVE** restituire il puntatore originale (`ExprPtr`).
- È **severamente vietato** effettuare "Deep Copy" o clonazioni ricorsive dell'albero se non strettamente necessario (es. cambio di Arena).
- L'identità del puntatore è il metodo primario per verificare l'uguaglianza in $O(1)$.

### 3. Vera Memory Arena (Bump Allocator)
- I nodi AST risiedono in una `AstArena` che alloca **blocchi contigui di memoria** (es. 64KB).
- È **vietato** usare `std::make_unique` o `new` per ogni singolo nodo (distrugge la cache L1/L2).

### 4. Moltiplicazione Implicita
- Il Parser **deve** iniettare implicitamente il token `*` tra Primary Expressions adiacenti (es: `2x` -> `2*x`).

### 5. Gestione DecimalLit (Confine Simbolico/Numerico)
- I letterali decimali (`DecimalLit`) sono ammessi solo per preservare l'input utente.
- Il Core Simbolico deve restituire `CASErrorKind::Unimplemented` se un'operazione algebrica (es. espansione) richiede di operare su un `DecimalLit`. Il calcolo numerico è un fallback esplicito.

---

## WORKFLOW: PRIMA DI DICHIARARE UN TASK COMPLETATO

1. **Integrità Matematica**: Tutti i test unitari e di integrazione passano al 100%.
2. **Benchmark Gate (OBBLIGATORIO)**: 
   - Esegui `bash scripts/benchmark.sh`. 
   - Se le prestazioni degradano rispetto alla baseline (`baseline_release.txt`), il codice va ottimizzato prima del merge.
3. **Sanitizers**: Zero errori in ASan e UBSan.
4. **Zero Warning**: Compilazione pulita con `-Wall -Wextra -Wpedantic -Werror`.
5. **Trace Validation**: Se il task tocca il `Simplifier`, verifica che il `ComputationTrace` sia accurato e non ridondante.

---

## STANDARD TECNICI E ANTI-DEBITO

- **Linguaggio**: C++20 rigoroso.
- **Anti-Monolito**: Limite massimo di **500 righe per file sorgente**. Oltre tale soglia, lo split in moduli specializzati è **obbligatorio** (es. `simplify_arithmetic.cpp`, `simplify_functions.cpp`).
- **Error Handling**: Divieto assoluto di `throw/catch` nel core. Solo `Result<T>` monadico. Ogni path di errore deve essere tracciabile via `CASError`.
- **Gestione Memoria**: 
  - Solo `AstArena`. Vietato l'uso di `std::shared_ptr` o `std::unique_ptr` per i nodi AST.
  - **Dangling Protection**: Vietato passare `ExprPtr` tra arene diverse a meno che non sia una copia esplicita (re-interning). Nei test, usa sempre l'arena del `CASContext` per il parsing.
- **Integrità Matematica**: 
  - Ogni regola di semplificazione deve essere orientata via **LPO (Lexicographic Path Ordering)** per prevenire cicli infiniti.
  - Niente "cerotti" (if-else hardcoded per casi base). Usa il sistema di rewrite universale.
- **Testing**: Solo confronto strutturale o equivalenza matematica. **Mai** usare `toString()` per validare la logica.
- **UI/UX Desktop**: Il tool `cas_ui` è un visualizzatore professionale basato su **ImGui**. Supporta campionamento adattivo e plotting 2D/3D (Fase 10).

---

## DOCUMENTAZIONE DI RIFERIMENTO
1. Regole Architetturali: `.APROJECT_REFERENCES/01_PROJECT_GOALS/02_architectural_rules.md`
2. Roadmap (FASE 1b): `.APROJECT_REFERENCES/02_ROADMAP/01_full_roadmap.md`
3. Plotting & Sampling: `.APROJECT_REFERENCES/03_ENGINE_MODULES/11_plotting_sampling.md`

---

## PROTOCOLLO DI AUTO-EVOLUZIONE (Self-Update)

L'IA ha il mandato di proporre aggiornamenti a questo file nei seguenti momenti chiave della timeline:

### 1. Transizione di Fase (Context Shift)
All'inizio di ogni nuova Fase della Roadmap (es. passaggio dalla Fase 1b alla Fase 2), l'IA deve:
- Aggiornare i link alla documentazione di riferimento nel paragrafo sopra.
- Aggiungere eventuali "Vincoli di Dominio" specifici della nuova fase (es. regole per la derivazione simbolica nella Fase 4).

### 2. Post-Mortem di un Bug Sistemico
Se durante lo sviluppo emerge un errore strutturale o un pattern che causa regressioni di performance:
- L'IA deve aggiungere una riga nella sezione "REGOLE ARCHITETTURALI" per proibire esplicitamente quel pattern (es. "Vietato usare X perché causa Y").

### 3. Ottimizzazione del Workflow
Se uno script di automazione (es. `benchmark.sh`) viene modificato o ne viene aggiunto uno nuovo:
- L'IA deve aggiornare la sezione "WORKFLOW" per includere il nuovo comando obbligatorio.

### Vincoli di Modifica
- **MAI** rimuovere o indebolire le Regole 1, 2 e 3 (BigInt, Structural Sharing, Arena). Sono le fondamenta immutabili.
- Ogni modifica al `CLAUDE.md` deve essere comunicata allo sviluppatore con la dicitura: *"Aggiornamento Costituzione Tecnica: [MOTIVAZIONE]"*.
