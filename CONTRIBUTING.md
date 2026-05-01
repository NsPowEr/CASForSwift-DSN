# CONTRIBUTIONS & CODING STANDARDS — REAL CAS ENGINE

Questo file sancisce la fine della Fase 0 e funge da monito per tutti i contributori e gli agenti IA.

## Le Red Lines dello Sviluppo (Mai Superarle)

1. **IL C++ È LA LEGGE FINALE**: Lo sviluppo si concentra unicamente sui sorgenti all'interno della cartella `src/` ed `include/`. NESSUNA UI in Swift deve essere creata in questo repository. 
2. **STANDALONE UI ISOLATA**: Tutti i test e le manipolazioni umane manuali devono avvenire lanciando ed eseguendo ESCLUSIVAMENTE la target locale allocata in `tools/cas_ui`.
3. **ZERO WARNINGS POLICY**: La build di CMake è stata configurata per rifiutare (`-Werror`) compilazioni difettose. Non patcheremo mai la policy, risolveremo sempre i warning.
4. **DOCUMENTAZIONE VANGELO**: La directory `.APROJECT_REFERENCES` istruisce al millimetro ogni agente. Nessun design o tool che non sia codificato e approvato lì dentro può approdare nel codebase.
5. **NESSUN MOCK**: I test matematici passano se e solo se la logica algebrica è oggettivamente e strutturalmente fondata (confronto degli alberi `ExprPtr`).

Leggere `CLAUDE.md` come entry-point globale.
