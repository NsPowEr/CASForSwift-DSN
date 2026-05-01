Piano di sviluppo per eliminare hardcode e migliorare il CAS engine

1 Contesto e obiettivi

Il progetto di un motore CAS per calcolatori simbolici presenta ancora lacune rispetto alla calcolatrice HP Prime G2.  La valutazione ACID evidenzia tre fallimenti principali:

1. Test 28 — limite \lim\limits_{x	o\infty} rac{\ln(\ln(x+e))}{\ln(\ln x)}: il motore usa la tecnica MRV (most rapidly varying) che sostituisce x=1/\omega e applica espansioni di Taylor. Per espressioni con logaritmi annidati la serie di Taylor in \omega ritorna erroneamente 0, impedendo la chiamata al metodo try_log_log_limit.
2. Test 9 — radici dell’equazione x^6-1: il codice attuale genera sei RootOf opachi. Il polinomio è un prodotto di fattori ciclotomici di grado ≤4;  l’algoritmo non tenta la fattorizzazione e non riconosce i fattori ciclotomici.
3. Test 19 — sistema di Gröbner: il file polynomial_groebner.cpp è ancora un stub; la funzione solve_nonlinear_system_f4 include un hardcode per un caso di test e ignora la base di Gröbner calcolata.

Oltre a questi punti, il piano REV 3 prevede l’implementazione di numerose funzionalità (ring di serie troncate, classificatore ODE, sommatoria ipergeometrica, algoritmo di Kovacic, ecc.) che sono ancora assenti.  L’obiettivo è definire un percorso di sviluppo approfondito e senza patch su misura: le funzioni dovranno utilizzare algoritmi generali ben documentati.

2 Analisi degli algoritmi di riferimento

2.1 Algoritmo F4 per basi di Gröbner

L’algoritmo F4 di Faugère è una variante efficiente dell’algoritmo di Buchberger.  Il principio base è trattare più coppie (S‑polinomi) contemporaneamente: si seleziona un sottoinsieme di coppie, si costruiscono polinomi ridotti e si inseriscono le combinazioni in matrici di Macaulay.  Queste matrici, che sono sparse, vengono poi ridotte in forma di Gauss (row‑echelon) per determinare nuovi generatori.  La descrizione di riferimento sottolinea due punti essenziali:

* F4 seleziona un sottoinsieme di coppie P’ e usa una funzione reduction che restituisce nuovi elementi di base.  La riduzione sostituisce la divisione polinomiale con la riduzione di righe di una matrice: i polinomi generati in symbolicPreprocessing vengono disposti come righe, la riduzione gaussiana produce combinazioni lineari ridotte, e i polinomi con monomi guida nuovi vengono aggiunti alla base .
* In symbolicPreprocessing vengono calcolati i contributi \mathrm{LCM}(\mathrm{LM}(g_i), \mathrm{LM}(g_j))/\mathrm{LT}(g_i) \cdot g_i e analoghi per g_j.  Quindi si aggiungono multipli di polinomi in modo che per ogni monomio divisibile da un monomio guida esista già la corrispondente combinazione nella matrice .  Questo assicura che le operazioni di divisione siano simulate da operazioni di riga.

Implementare F4 significa separare chiaramente:

1. Selezione e filtraggio delle coppie con criteri di produttività (ad es. evitare S‑polinomi che si annullano per il product criterion).
2. Preparazione simbolica: generare le righe della matrice Macaulay includendo tutte le combinazioni necessarie.
3. Riduzione di matrice: usare eliminazione gaussiana o fattorizzazione LU sui sistemi sparsi; l’efficienza dipende dal trattamento di matrici grandi ma sparse .
4. Inter‑riduzione: una volta ottenuta la base, ridurre ciascun polinomio rispetto agli altri per ottenere una base ridotta (unica).

2.2 Polinomi ciclotomici

Per n≥1, il polinomio ciclotomico \Phi_n(x) è l’unico polinomio monico con coefficienti interi che divide x^n-1 ma non divide x^k-1 per k<n.  I suoi zeri sono le radici primitive n‑esime dell’unità  e soddisfa la relazione
\prod_{d\mid n}\Phi_d(x) = x^n-1.
Il teorema di Möbius fornisce una formula esplicita per \Phi_n(x):
\Phi_n(x)=\prod_{d\mid n}(x^d-1)^{\mu(n/d)},
con \mu funzione di Möbius .  Questa formula consente di calcolare \Phi_n ricorsivamente dividendo x^n-1 per i \Phi_d con d<n  e permette di costruire una tabella delle prime 30 \Phi_n per riconoscere fattori ciclotomici di grado piccolo.

2.3 Limiti simbolici e serie gerarchiche

L’algoritmo MRV/Gruntz valuta i limiti trasformando le espressioni in forma standard con x→∞, identifica l’insieme delle sottoespressioni più rapidamente varianti tramite il limite del rapporto dei logaritmi  e usa una variabile piccola \omega = 1/x per espandere la funzione in serie: i termini indipendenti da \omega determinano il limite .  Tuttavia, quando compaiono logaritmi annidati come \ln(\ln(x+e)), la sostituzione x=1/\omega crea una singolarità essenziale in \omega=0; l’espansione di Taylor restituisce 0 e l’algoritmo termina prematuramente (come osservato nel Test 28).

Un’alternativa più robusta è l’uso di serie Puiseux/hierarchiche, che permettono esponenti frazionari e coefficienti contenenti logaritmi annidati.  Nella letteratura sui sistemi algebrici si sottolinea che tali serie gestiscono correttamente singolarità essenziali e logaritmi annidati : gli algoritmi di Geddes, Gonnet e Gruntz producono serie troncate in cui i termini sono ordinati per “dominanza” e rispettano le dipendenze gerarchiche fra logaritmi e potenze.

2.4 Equazioni differenziali lineari a coefficienti costanti

Per un’equazione omogenea
y'' + p\,y' + q\,y = 0,
la soluzione si ricava risolvendo l’equazione caratteristica r^2 + pr + q =0.  A seconda del discriminante \Delta = p^2 - 4q, si ottengono tre casi :

* \Delta>0: due radici reali r_1,r_2 ⇒ y = C_1\,e^{r_1 x} + C_2\,e^{r_2 x}.
* \Delta=0: radice doppia r ⇒ y = (C_1 + C_2 x) e^{r x}.
* \Delta<0: radici coniugate v\pm w\,i ⇒ y = e^{v x}\,igl(C_1\cos(w x) + C_2\sin(w x)igr).

Per un’equazione non omogenea y'' + p y' + q y = f(x) il metodo della variazione dei parametri fornisce una soluzione particolare:
y_p = -y_1\int rac{y_2\,f}{W}\,dx + y_2\int rac{y_1\,f}{W}\,dx,
con y_1,y_2 soluzioni della parte omogenea e W=y_1 y_2'-y_2 y_1' il Wronskiano .

2.5 Algoritmo di Gosper per sommatorie ipergeometriche

Data una sequenza ipergeometrica s_k (cioè con rac{s_{k+1}}{s_k} razionale), Gosper cerca una sequenza z_k tale che s_k = z_{k+1} - z_k.  La lezione di algebra computazionale mostra che il problema si riduce a trovare una funzione razionale y_k con z_k = y_k s_k e risolvere una ricorrenza linearea_k y_{k+1} - y_k = 1 dove a_k=rac{s_{k+1}}{s_k} .  La ricorrenza può essere ulteriormente semplificata scrivendo a_k = rac{p_{k+1} q_{k+1}}{p_k r_{k+1}} con p,q,r polinomi a coprimezza controllata, e introducendo un polinomio f_k tale che z_k = (r_k/p_k) s_k f_{k-1}.  Sostituendo in s_k = z_{k+1} - z_k si ottiene un’equazione lineare con coefficienti polinomiali per f_k ; risolvere questo sistema fornisce l’antidifferenza z_k e quindi la somma telescopica.

3 Pianificazione delle modifiche

3.1 Fase 2b: limiti con logaritmi annidati (Test 28)

1. Chiamata anticipata a try_log_log_limit.  Prima di invocare l’algoritmo MRV, verificare se l’espressione ha la forma \ln(A)/\ln(B) e se \lim_{x	o\infty}rac{A}{B} è finito e non nullo.  In tal caso il limite è 1.  Occorre spostare questa logica prima della sostituzione x=1/ω per evitare la serie di Taylor difettosa (il piano REV 3 suggerisce di inserire un blocco di sei righe nel file limit.cpp).
2. Estensione di compare_growth.  Attualmente compare_growth conta solo la profondità di exp().  Per logaritmi annidati è necessario definire una metrica che consideri anche la profondità di ln().  Una possibile estensione consiste nel definire classi di comparabilità basate su \lim_{x	o\infty} rac{\log|f|}{\log|g|} (come nella definizione MRV ) ma applicata ricorsivamente a funzioni logaritmiche.  Se l’espressione contiene \ln(\ln(x)) e \ln(\ln(x+a)), tali componenti devono comparire nella stessa classe e la sostituzione x=1/ω non deve essere applicata.
3. Fallback con serie gerarchiche.  Nel caso in cui l’algoritmo MRV non riesca a determinare il limite (ad esempio quando l’espansione su ω produce 0), implementare un modulo series_expand che costruisce una serie Puiseux o gerarchica.  Le serie gerarchiche permettono coefficienti contenenti logaritmi annidati e ordinano i termini in base alla dominanza .  Questo modulo dovrebbe supportare operazioni di base (somma, prodotto, composizione, inversione) e funzioni elementari (exp, ln, potenze).  Integrare series_expand in limit_series.cpp come fallback quando l’approccio derivativo fallisce.
4. Test di regressione.  Dopo l’implementazione eseguire i limiti di riferimento (Gruntz su x→0, Schanuel, Squeeze) per assicurarsi che l’anticipo della regola log/lim non introduca regressioni.

3.2 Fase 9: radici ciclotomiche (Test 9)

1. Estendere solve_polynomial.  Nella funzione di soluzione dei polinomi, prima di emettere RootOf per gradi ≥5, tentare sempre la fattorizzazione completa mediante solve_by_factoring.  Nel caso x^6 - 1 la fattorizzazione produce quattro fattori di grado ≤2; le soluzioni si ottengono risolvendo questi fattori.
2. Riconoscimento ciclotomico.  Implementare un modulo polynomial_cyclotomic.cpp che costruisca \Phi_n(x) per n≤30 con la formula di Möbius .  Fornire una funzione is_cyclotomic(const IntPoly&) → optional<int> che, dato un polinomio intero monico, restituisce n se p=\Phi_n.  Usare questa funzione in solve_factor per riconoscere fattori ciclotomici; quando p=\Phi_n, emettere soluzioni esplicite e^{2\pi i k/n} per k coprimo con n.
3. Integrazione con la fattorizzazione.  solve_by_factoring deve invocare factorization_polynomials per decomporre il polinomio in fattori primi e poi chiamare ricorsivamente solve_factor su ciascun fattore.  Se un fattore è ciclotomico, generare soluzioni; altrimenti, se il grado ≤4, usare le formule generali.  Solo se la fattorizzazione fallisce, usare RootOf.
4. Test e ottimizzazioni.  Verificare le seguenti istanze: solve(x^6-1,x) produce sei radici evidenti, solve(x^4-1,x) produce {1,−1,i,−i} e solve(x^2+x+1,x) produce le due radici complesse di \Phi_3.  Per gradi più alti, assicurarsi che la ricerca ciclotomica non rallenti eccessivamente il percorso.

3.3 Fase 1: integrazione completa dell’algoritmo F4

Il file polynomial_groebner_f4.cpp contiene già un’implementazione parziale di F4 (dati i tipi PolyF4 e MacaulayMatrix).  Restano da completare diversi punti.

1. Separazione in header.  Creare polynomial_groebner_f4.hpp che dichiari le strutture (Monomial, PolyF4, Pair) e le funzioni f4_groebner, f4_to_expr, expr_to_f4.  Definire MonomialOrder con le varianti Lex e GRevLex e il comparatore MonomialGRevLexComparator.  Questo header deve essere incluso sia in polynomial_groebner_f4.cpp sia in polynomial_groebner.cpp.
2. Parametrizzazione sull’ordine.  Modificare f4_groebner affinché accetti un parametro MonomialOrder order.  In base all’ordine, istanziare MacaulayMatrix con il comparatore appropriato.  Implementare il product criterion di Buchberger per evitare coppie con monomi guida coprimi.
3. Conversione da AST a PolyF4.  Implementare expr_to_f4(ExprPtr, vector<Symbol>, CASContext&) che costruisca PolyF4 direttamente dall’AST usando coefficienti razionali.  L’algoritmo ricorsivo dovrà gestire costanti, razionali, simboli, somme, prodotti, potenze intere e n‑arie.  Questa funzione permetterà di bypassare le limitazioni di MultivariatePolynomial e di accettare coefficienti frazionari.
4. Inter‑riduzione.  Aggiungere una routine inter_reduce(vector<PolyF4>&, MonomialOrder) che riduca ogni polinomio della base rispetto agli altri: se esiste un monomio in G[i] divisibile dal monomio guida di G[j], sottrarre la combinazione adeguata e rendere G[i] monico.  La base ridotta garantisce determinismo.
5. Wrapper polynomial_groebner.  Riscrivere completamente polynomial_groebner.cpp facendo: (1) convertire le equazioni in PolyF4 tramite expr_to_f4, (2) chiamare f4_groebner con GRevLex per efficienza, (3) inter‑ridurre la base, (4) convertire i polinomi risultanti in ExprPtr con f4_to_expr.  La funzione deve restituire l’intera base di Gröbner in Lex o GRevLex, a seconda dell’overload.
6. Soluzione di sistemi polinomiali.  Per solve_nonlinear_system_f4 eliminare il codice hardcoded.  Implementare il shape lemma: calcolare la base di Gröbner in ordine Lessicografico, ridurre la base, individuare il polinomio “puro” in una singola variabile e risolverlo tramite solve_polynomial.  Effettuare back‑substitution ricorsiva per ottenere tutte le soluzioni.  Gestire casi di dimensione zero (ideale 0‑dimensionale) e restituire “Unimplemented” altrimenti.

3.4 Fase 4: soluzione di sistemi non lineari via Gröbner

Una volta completata l’implementazione F4, modificare csolve.cpp affinché utilizzi la base di Gröbner per risolvere sistemi n×n:

1. Se il sistema è 2×2, mantenere il percorso esistente via risultato (subresultant).
2. Altrimenti, chiamare polynomial_groebner con ordine Lex per ottenere una base di Gröbner.  Applicare lo stesso algoritmo di back‑substitution del solve_nonlinear_system_f4.  Restituire le soluzioni come array di valori.

3.5 Fase 6: classificazione ODE e risolutore a coefficienti costanti

1. Classificatore.  In ode_classifier.cpp implementare la ricerca di equazioni del tipo a2(x) y'' + a1(x) y' + a0(x) y = f(x) spostando tutti i termini a sinistra e raccogliendo i coefficienti mediante collect_coefficient.  Usare is_constant_wrt(expr,x) per verificare se a2,a1,a0 sono costanti; se sì, classificare come Linear2ndOrderConstantCoeff e fornire i componenti [a2,a1,a0,f].  Altrimenti tentare la classificazione di primo ordine (lineare, Bernoulli, separabile) e restituire Unknown altrimenti.
2. Risolutore avanzato.  In ode_solver_advanced.cpp, dato un classificatore Linear2ndOrderConstantCoeff, calcolare il discriminante \Delta.  Usare i casi descritti sopra  per costruire y_h.  Per il termine particolare, implementare il metodo della variazione dei parametri: calcolare il Wronskiano W e le integrali \int (\pm y_i f)/(a_2 W)\,dx .  Se l’integrale non è risolvibile, restituire solo la soluzione omogenea con un avviso.
3. Test.  Eseguire le equazioni campione: y'' - y = 0, y'' + 4y = 0, y'' - 3y' + 2y = 0, y'' + y = \sin x, y'' + 2y' + y = 0 e confrontare con le soluzioni note.

3.6 Fase 2c: serie troncate e gerarchiche

1. Struttura dati.  Implementare TruncatedSeries come descritto nel piano REV 3: variabile, centro, ordine di troncamento, valutazione (esponente minimo), denominatore Puiseux, vettore di coefficienti.  Consentire esponenti negativi e frazionari.
2. Operazioni del ring.  Implementare funzioni per somma, prodotto, negazione, inversione (Newton), composizione e derivazione.  Per funzioni trascendenti (exp, log, potenze, sin, cos) implementare algoritmi basati su ODE o tabelle di Taylor.  Assicurarsi che la serie logaritmica richieda a_0=1 e l’esponenziale a_0=0.
3. Espansione da AST.  Scrivere series_expand(expr,t,center,order,ctx) che analizza l’AST e produce una TruncatedSeries usando le operazioni sopra.  Integrare series_expand in limit_series.cpp e taylor_series come fallback quando la derivazione ripetuta fallisce.
4. Test: verificare le espansioni 1/(1-x) (serie geometrica), 	an x, \exp(\exp(x)-1); confrontare con i risultati del pianificatore HP Prime.

3.7 Fase 7: sommatoria simbolica (Gosper e Zeilberger)

1. Riconoscimento ipergeometrico.  Implementare hypergeometric_ratio(t_k,k,ctx) che calcola r(k)=t(k+1)/t(k) razionalizzando la frazione; restituire numeratore e denominatore come polinomi in k.  Se la funzione non è razionale (ad es. t_k=k!\cdot \sin k), restituire un errore.
2. Algoritmo di Gosper.  Sviluppare gosper_antidifference(t_k,k,ctx) seguendo lo schema descritto: (1) calcolare il fattore ratio r(k)=p(k+1)q(k+1)/(p(k) r(k+1)); (2) trovare i polinomi a,b,c (fattorizzazione di Gosper) con condizioni di coprimezza ; (3) stabilire un limite superiore sul grado del polinomio f(k) e risolvere l’equazione funzionale a(k)f(k+1)-b(k-1)f(k)=c(k) in f(k) ; (4) ricavare z_k e quindi la somma.  Se il sistema lineare non ha soluzione, restituire Unimplemented.
3. Algoritmo di Zeilberger.  Per somme con parametro, implementare una versione semplificata: per una funzione ipergeometrica F(n,k) verificare la razionalità di F(n,k+1)/F(n,k), poi cercare una ricorrenza di ordine J ≤3 nella variabile n usando Gosper parametrico.  Per casi semplici (somme polinomiali, geometriche), utilizzare formule chiuse (polinomi di Bernoulli, formule geometriche).  Per somme definite, usare il teorema fondamentale del calcolo discreto: \sum_{k=a}^b t_k = Z(b+1) - Z(a).
4. Integrazione in summation.cpp.  Sostituire il blocco Unimplemented con: tentativo di Gosper per somme indefinite, poi per somme definite, poi riconoscimento polinomiale, infine Zeilberger per somme con parametri.

3.8 Fase 6b: algoritmo di Kovacic (L=1)

1. Analisi dei poli di r(x).  Per l’equazione y'' = r(x) y con r \in \mathbb{Q}(x) determinare i poli di r.  Per ogni polo c di ordine n calcolare le possibili parti principali di \omega(x), soluzione dell’equazione di Riccati \omega'+\omega^2=r.
2. Costruzione di \omega.  Scrivere \omega(x) = A(x) + \sum 	ext{parti principali} dove A(x) è un polinomio di grado massimo d (calcolato in base ai poli).  Sostituire in \omega'+\omega^2=r e risolvere per i coefficienti di A con un sistema lineare (usando l’algebra lineare esistente).
3. Soluzione.  Se la soluzione esiste, integrare \omega e costruire y = \exp(\int \omega\,dx).  Se nessuna soluzione viene trovata, restituire Unimplemented.  Assicurarsi che l’algoritmo chiami i moduli di fattorizzazione e integrazione razionale già implementati.

3.9 Infrastruttura e test

* Riduzione di Unimplemented: contare le occorrenze di Unimplemented prima e dopo ogni fase per assicurarsi che diminuiscano, come indicato negli script verify_milestone.sh.
* Suite di test: creare test unitari dedicati (GB, ODE, sommation, serie) e aggiornare test_hpprime_acid.cpp per includere i nuovi casi.  Assicurarsi che i tempi di esecuzione rimangano sotto i limiti dei benchmark.
* Benchmarking: eseguire scripts/benchmark.sh dopo ogni fase e confrontare con baseline_release.txt.  Ottimizzare eventuali regressioni di performance.

4 Conclusioni

Il programma di lavoro delineato richiede un impegno significativo ma porta a una struttura robusta e modulare.  Ogni funzione nuova o corretta deve basarsi su algoritmi noti: la riduzione di polinomi mediante F4 e inter‑riduzione, il riconoscimento dei fattori ciclotomici con la formula di Möbius, l’espansione in serie gerarchica per limiti difficili, la classificazione formale delle ODE con soluzioni standard, e la sommatoria tramite algoritmi di Gosper e Zeilberger.  Seguendo questa roadmap, il motore CAS potrà eliminare i hardcode e avvicinarsi alla parità funzionale con la HP Prime G2.