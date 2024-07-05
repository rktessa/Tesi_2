# Master Thesis repository
## Understanding the effect of opinions and behaviours on the spread of infectious diseases
### Riccardo Tessarin 2023/2024

Development of a multilayer system model capable to explain some functional shape recurring in behaviour dynamics influenced by an epidemic spread.

Two system layer:
 - SIRS for epidemic;
 - Careless, Compliant, Against (double contagion from initial unaware population) behaviour model.

link for data:

https://ourworldindata.org/covid-income-support-debt-relief

Note incontro:

- Nuova struttura di divisione della tesi, ipotesi di spostare il capitolo 2 e 3 unendoli e lasciare il 4 come "development" e non "description". Dividere la parte in cui presento teoria e lavori fatti da altri, a quella in cui è presentato quello che sto sviluppando.
- Andare avanti con lo studio di analisi del modello doppio. una idea potrebbe essere quella di cercare di ridurre la dimensionalità del sistema, e fissare i parametri che possono considerarsi più stabili, come ad esempio la fatigue,che sembra dai dati abbastanza piccola. Il libro che consigliava Daniele, contiene delle sezioni in cui vengono affrontate queste tecniche. Leggi e copia, fai, non leggere e basta.
- rispetto anche a lavoro simile nel modello come quello dell'articolo, si vede una certa differenza di approccio, obbiettivi e prospettiva. Si può spiegare anche quello e tenere come reference metodologica.
- il mio lavoro ha come obbiettivo avere una attinenza a dei dati realmente osservabili, 
- dai dati vedere come riproduco i pattern e stima dei valori dei parametri
- analisi di sensitività serve per trovare delle stime in cui hanno un effetto i vari parametri e poi quelli saranno i boundaries da cui eseguire analisi di fitting.
- io sto descrivendo due modelli piccoli e poi li unisco
- sempre per come ci si pone con la letteratura capire che se anche c'è già qualcosa di fatto, ma il mio approccio è nuovo è quello che conta. ed è quello che sto facendo ora.  
- bisogna rivedere codice del modello misto e poi correggere il delta con quello che è corretto e cambiare il nome di uno dei due
- all'inizio, inquadrare la situazione, cioè le epidemie esistono, poi alllora esprimere il desiderio di studiare questo fenomeno e come venga fatto bene da tempo. che il modelling che volgio fare aggiunge anche degli altri obiettivi multidisciplinari, volendo integrare l'effetto del comportamento delle persone. e l'awareness  che è un filone molto attivo e promettende nowadays. e può essere utile per modelli che permettenao di capire meglio le epidemie e pe queloso la tesi si vuole inserire in questo filone di ricerca. 

- dopo la parte di descrizione dei vari modelli inserire un capitolo di objectives e di summary dei contenuti. nella prima parte ripeterai un pochino quello che hai detto nella introduction e poi descirvi i vari capitoli. 
- per lo studio del modello totale fai riferimento anche a come è stato analizzato quello della onchcerchiasis. 
- l fatigue varia poco dallo scraping di daniele, da individuare i parametri più sensibili dei due modelli separati. e provare poi ad estendere, nel modello complessivo e sperare che le cose correspondano come effetto. 



## Cose da finire il 9 maggio  
Gli articoli che devo finire di guardare per la parte sir sono: 
la review: per gli agent based
appunti da daniele: per i multi layer syste/network

articolo hetcote se serve ancora
articolo okabe per altre info utili su sir

## Significato del modello
Dividere la popolazione in 
- Compliant/Carefull
- Heedless/Careless
- Against rules/Risky

A un tempo zero la popolazione è SH nella quasi totalità. Andando avanti si dividerà in SC e SA e poi questi tratti rimangono negli strati sottostanti, perchè a quel punto più il tempo va avanti più le persone hanno una memoria storica e ha senso che rimangano i due gruppi pro e contro. Quindi il modello tiene conto a questo livello di S dell'esistenza di una zona grigia in cui stanno le persone senza una opinione forte. A livello di una situazione in cui la pandemia è appena iniziata mi aspetto siano tutti qui e anche che sia un comparto che si riempie abbastanza quando il disease c'è  da più tempo, continua a girare ma la tipo è meno mortale. All'inizio non c'è timore, poi si diffonde una ‘‘concerned awareness,’’!

COme diceva la prof le person IC sono più a casa e quindi partecipano meno all'infezione (anche sensibilmente meno forse, poi questi coefficienti vanno tarati).


Nella accezione di behaviour avevo una idea rispetto al senso... 
- i compliant stanno a casa e non partecipano all'nfezione mentre gli against non rispettano queste regole.

- i compliant si proteggono, gli against no.
 - Decidere che modello usare per awareness della popolazione. 

 In epi_behaviour01.m ho tolto la mortalità natalità come effetto, mentre ho cambiato i nomi delle variabili e corretto il codice per adeguarlo alle decisione nuove prese.


21 Giugno:

- Vedere dove eri arrivato a scrivere prima per concludere beneddetto primo capitolo. 
- Sempre nel primo capitolo hai trovato un po' di materiale da inserire per la parte di glossario di Opinion e Behavior. 

### Feedback Daniele 

 Ho alcune domande e suggerimenti. Inizio a scriverteli qui, poi magari ci si può sentire a inizio settimana prossima per approfondire:

    1. Phi è di fatto un parametro moltiplicativo, che tu interpreti come awareness (ma che forse è più direttamente legato alla "forza" di attrazione di un gruppo, è una sorta di parametro moltiplicativo di mixing). Inoltre, che questa "awareness" o "attrattività" sia uguale a prescindere dallo stato lo vedo un'ipotesi audace: potrebbe benissimo essere state-dependent (es: un infetto può "pentirsi" e diventare compliant più facilmente). Purtroppo non ci sono studi psicologici sul campo - o non sono riuscito a trovarne...
    2. In buona sostanza: forse attenuerei interpretazioni così forti e riassorbirei Phi in k_3, che rimane così un parametro di contatto tra categorie che la "pensano" diversamente e si influenzano a vicenda. Di fatto, i tuoi risultati non cambiano di una virgola (si tratta solo di cambiare lettera), ma non si rischia di scivolare su interpretazioni ardite. Viceversa, i k_i saranno interessanti per valutare la forza e l'impatto della peer pressure (magari vs un "controllo mean field" di cui avevamo parlato).
    3. Nei layer I e R ci sono flussi orizzontali dovuti sia a "fatigue" che a peer pressure. Una buona idea: aumenta i parametri, ma consente di fare studi sul ruolo della peer pressure stessa :)
    4. A proposito di interpretazioni, non sono convintissimo di come hai spiegato Phi_n: così com'è ora è una sorta di probabilità, ma pensiamo bene a come chiamarlo
    5. In (3.5), non dovresti usare x,y,z definiti prima, per i compartimenti? Oppure ridefinisci le variabili (oppure c'è un passaggio che manca)
    6. \beta nelle equazioni per la next gen matrix: usa un altro simbolo (anche solo \beta') per non confonderti con quello del modello
    7. (3.6): é il risultato esatto o dovresti anche invocare il leading eigenvalue del prodotto matriciale? Dipende da come ha definito le equazioni; non ho guardato in dettaglio, prova a ridare un'occhiata e sii il più generico possibile (o spiega perchè va bene così nel caso specifico)
    R0: scrivilo come formula e non come immagine. Non si riesce a compattare meglio? Magari già usando le considerazioni sui valori di A0 e C0 che usi nel paragrafo sotto?
    8. Case I: deriva la formula generica prima di sostituire valori specifici. Puoi anche fare considerazioni di scala e sostituire alcuni valori -> 0 (nella parte analitica). Riesci a tirare fuori dei range di valori di determinati parametri, tali per cui R0 < 1 sempre? E a calcolare analiticamente la sensibilità ad alcuni parametri, usando misure di sensitivity analysis quali |\partial R0 / \partial p_i| (dove p_i è ogni parametro)?
    9. Ci sono poi alcune cose da sistemare nel testo e nella spiegazione per evitare di lasciare punti aperti o consecutio non chiarissime, ma in generale direi che il discorso fila bene! Per quando ci sentiamo, prova a buttare giù una lista di investigazioni che ti piacerebbe fare per proseguire, e ne parliamo.
    Riesci a fare un disegno schematico del modello? Aiuterebbe a capirlo meglio 🙂


Conti:
- cercato di fare uno scaling e arrivato a qualcosa. Capire come andare oltre. tolto una equzione e riscritto in fuzione di meno parametri per ora. 

- voglio riscrivere in bella i conti per R-0-epi-beha e trovare una forma che toglie beta/gamma così plotto in funzione dei parametri di comportamento (e di gamma) e osservo cosa è più importante lì.




