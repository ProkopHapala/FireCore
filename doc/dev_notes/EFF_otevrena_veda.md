

# Projekt Otevrena Veda

# Úvod / Pozadní

* Elektronove silove pole ktere nahrazuje kvantove interakce mezi elektrony jednoduchym modelovym potencialem mezi (1) elektronem-a-elektronem a (2) elektronem a jadrem, a to za zjednodušeného předpokladu že vlnová funkce každého elektronu $i$ je gaussovský obláček $\psi_i(r)=\sqrt{C}\exp(-(|r-R_i|/w)^2)$ a jeho hustota pravdepodobnosti vyskytu $\rho_i(r)=C\exp(-2(|r-R_i|/w_i)^2)$ (kde C je normalizacni konstanta tak aby celkova pravdepodnost  byla 1 $\int \rho_i(r)=1$) 
   * Silové pole má oproti kvantovým metodám obroskou výhodu v rychlosti a názornosti.
   * Vetsi cast tohoto siloveho pole je odvozena rigorozně z prvních principů (Kinetická energie a Elektrostatická interakce)
   * Co není odvozeno rigorozně je tzv. *Pauliho repulze* mezi elektrony která je způsobena tím že eketrony jsou fermiony a dva fermiony nemohou okupovat stejný kvantový stat. 
     * Pzn. rigorozně vzato to vlastně není repulze, ale jen důsledek faktu že elektrony jsou nerozlišitelné, a pokud by byly ve stejném kvantovém stavu tak se jedná o elektron jeden (jakoby splynuly). Ale protože my chceme popisovat systém klasicky, s rozlišitelnými elektrony, tak se to projevuje jako repulze mezi nimi.
     * Tato energie se dá odhadnout z ortogonalizační energie mezi vlnovými funkcemi dvou elektronů, ale tato energie silně závisí na tvaru takže předpoklad kde elektronům vnutíme tvar gaussových orbitálu je příliž velkou aproximací, a takto odvozená Pauliho repulze neodpovídá úplně realitě (velká chyba).
     * Takto energie hraje roli také pro integraci valenčních elektronu s jádrem (resp. s iontem tvořeným jádrem stíněným vnitřními slupkami), protože zde hraje roli interakce s elektrony vnitřních slupek. 

# Cíl / myšlenka projektu

 * Cílem je najít vhodnou opravu k funkci popisující interakci valenčních elektronu se stíněným jádrem a dvou elektronů navzájem, tak aby výsledky co nejpřesněji popisovaly skuteční chování molekul - to znamená, aby molekuly relaxovaly do správné geometrie a aby závilost energie a síly na poloze atomů (Potential energy surface, PES) odpovídala přesnějším kvantovým metodám (density functional theory, DFT). 
 * Toho dosáhneme tak že vyjádříme opravu k funckím $$E_{ae}(r_{ae},s_e)$$ a $$E_{ij}(r_{ij},s_i,s_j)$$ v nějaké vhodné reprezentaci a budeme 
(např. součet Gaussových funkcí, polynom, Reciprokou funkci, Spline) popsanou nějakými parametry (např. rozvojové koeficient, poloměry gaussianů atd.) a budeme hledat takové parametry $\vec p$, které minimalizují rozdíl energie $\min_{\vec p}\sum_k(E_{EFF}(\vec r_k,\vec p) - E_{DFT}(\vec r_k)$ pro všechny trenovací geometrie molekul $r_k$  předpovězené takto opraveným Elektron Force Fieldem $E_{EFF}$ a refereční energie získané z DFT $E_{DFT}$.

# Problém / Výzva

* Hlavním problémem projektu je to že energie v elektron forcefield závisí ne jen na polohách jader ale také na polohách $R_i$ a poloměrech $w_i$ jednotlivých elektronů. Ty ale nelze přímo porovnat s kvantovými výpočty, protože v nich jsou elektrony delokalizovány (typicky v několika atomech a vazbách - viz. Molekulární orbitaly). 
   * => S $R_i, w_i$ tedy musím pracovat jako se skrytými parametry, což komplikuje použít tradiční metody fitování (řetězové pravidlo pro derivace, ze kterého byla odvozena back-propagace v machine-learningu). Misto toto budeme prinejmensim ze zacatku pouzovat Monte-Carlo fitovani kde budeme paramtery menit nahodne v okoli zatim nalezeneho nejlepsiho reseni (to je sice mene efektivni, ale zase nam to usetri praci a komplicace s vypocty derivaci energie podle parametru).
   * Jenou z možností je použít dodatečné regularizační vazby které budou penalizovat parametry opravy forcefieldu které vedou k hodnotám $R_i, w_i$ neodpovídajícími chemické intuici. 
     1. Ze základní organické chemie například víme že existuji $sp$,$sp^2$,$sp^3$ konfogurace atomů které mají lineární, trojuhelníkovou a tetragonální geometrii. Můžeme tedy měřit odchylku polohy elektronů a od těchto předpokládaných směrů a postihovat te které ji mají vyžší.  
     2. Poloměr jednotlivých elektronu můžeme zase odhadnout z dosahu elektronové hustoty okolo molekul, kterou dokážeme vypočítat z DFT.
   * Rigoroznější metoda by bylo přímo fitovat elektronvou hustotu z DFT pomocí rozvoje elektronového forcefieldu  $\min_{\vec p}\sum_k( \sum_i \rho_i(\vec r_k,\vec p) - \rho_{DFT}(\vec r_k))$
* Druhy problem / Otazka je v jake forme reprezentovat

# Postup / Kroky

* Vygenerujeme velke mnozstvi konfiguraci malych molekul $\vec r_k$ a spocitame je pomoci DFT
   * Ulozime: 
      * Energie
      * Sily na atomy
      * Atomove naboje (dipoly?)
      * Elektronovou hustotu - (v nejake vhodne reprezentaci - na mrizce, ve vyznecnych bodech ?)
          * Volba vhodne reprezentace hustoty muze byt jeden sub-problem
   * Pouzijeme asi FireBall (je mene presny ale dokazeme jeho kod upravovat jak potrebujeme). Pozdeji muzeme zpresnit napr pomoci psi4.
* Vybereme vhodnou reprezetnaci funkci 
   * $$E_{ae}(r_{ae},s_e)$$ - interakce elektronu s jadrem , dulezitejsi protoze urcije typ atomu
   * $$E_{ij}(r_{ij},s_i,s_j)$$ - interakce dvou elektronu navzajem - univerzalni, ale narocnejsi
* Zrelaxujeme $R_i, w_i$ elektronu pro kazdou trenovaci geometrii $\vec r_k$ Vykreslime si histogram jak nami generovane geometrie vzorkuji prostor.
   * Z toho zjistime co jsou relevanti geometrie, ktere mista na funkcich $$E_{ae}(r_{ae},s_e)$$, jestli je musime vzorkovat husteji
    
# Kratkodoby vyhled

* Zkusit si vykreslit funkci $$E_{ae}(r_{ae},s_e)$$ v prostoru - je to funkce 2 parametru ($$E_{ij}(r_{ij},s_i,s_j)$$ je funkce 3 parametru to je narcnejsi)
* napsat si jednodduchy Monte-Carlo optimalizacni algoritmus v pythonu



     
# English Background

Electron Force-field fitting consist of two parts

1. Electron-Core interaction $$E_{ae}(r_{ae},s_e)$$ where $s_e$ is size of the electron blob (gaussian with), and $r_{ae}$ is the distance between electron and the nuclei.

2. Electron-Electron interaction $$E_{ij}(r_{ij},s_i,s_j)$$ depending on sizes $s_i,s_j$ of the two electrons $i,j$ and their distance $r_{ij}$. 

* From both of these function we should subtract uncontroversial therms like coulomb interaction which is physically defined by error-function $E^{COUL}(r_{ij},s_i,s_j) = erf(r,s_{ij})= \int_r \exp(-(r_{ij}/s_{ij})^2)$ where $s^2_{ij}=s_i^2+s_j^2$.
   * $\Delta E_{ae}=E_{ae} - E^{COUL}_{ae}$
      * $=E_{ae} - erf(r_{ae}, s_{e})$
   * $\Delta E_{ij}=E_{ij}-E^{COUL}_ij$ 
     * $=E_{ij}-erf(r_{ij},s_{ij})$
* These are two scalar functions depending on 2-3 scalar arguments, therefore these function can be easily tabulated and plotted in 2D resp. 3D.
* While Electron-Core interaction $\Delta E_{ae}$ is simpler (just 2D) it depends on the type of atoms (charge of nucleus or ion). Generally this function have a role of **pseudo-potential** as it incorporates interaction of valence electron with both nuclei and core electron shells.
   * This makes $\Delta E_{ae}$ the key quantity if we want to tune EFF to reproduce Poential energy surface (PES) of different molecules, since we can choose different $\Delta E_{ij}$ depending not only on element (H,C,O,N) but also on refined atomic types (e.g. H_COOH, N_pyridine, etc.)
* The Electron-Electron interaction $\Delta E_{ij}$ should be universal, in principle. Also tuning this larger 3D function would be more challanging (as it requires more parameters to describe and more training examples to cover). 
   * However, if needed it would be possible to modify this function slightly depending on type of nearby atoms.
      * For example, assuming this interactin comprise mainly of Pauli exclusion and exchange-interactoin, we can estimate the it is dominated by the center-point between the two electrons (as densities of both gaussian blobs decay rapidly). 
      * physically we can argue that modification of $\Delta E_{ij}$ is due to interference with other electrons at the location. The property of electron cloud at the location can be approximated by properties of the atom.
      * We can then evaluate correction ot this function from interpolating nearby atomic nuclei. Nevertheless such operatin is definitely much more costly, than incorpoeratin of atomic types in $\Delta E_{ae}$.   