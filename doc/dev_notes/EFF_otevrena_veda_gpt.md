# Projekt Otevřená Věda

## Úvod / Pozadí

Elektronové silové pole (Electron Force Field, EFF) je model, který nahrazuje pravé kvantové interakce mezi elektrony jednoduchým potenciálem mezi gaussovskymi oblaky elektronové hustoty (které jsou mnohem jednodusí, než skutečné vlnové funkce elektrony). Gaussovská vlnová funkce každého elektronu je definovan jako:

$$ψ_i(r) = \sqrt{C} \exp\left(-\left(\frac{|r - R_i|}{s_i}\right)^2\right),$$

a hustota pravděpodobnosti jeho výskytu:

$$ρ_i(r) = |\psi_i(r)|^2 = C \exp\left(-2 \left(\frac{|r - R_i|}{s_i}\right)^2\right),$$

kde $R_i$ jsou polohy a $s_i$ poloměry elektronů, $C$ je normalizační konstanta tak, aby celková pravděpodobnost byla 1 ($\int ρ_i(r) \, dr = 1$).

 Tento potenciál popisuje:

1. Interakci mezi elektronem a jádrem $E_{ae}(r_{ae}, s_e)$
2. Interakci mezi dvěma elektrony $E_{ij}(r_{ij}, s_i, s_j)$


### Výhody a nevýhody EFF

- **Výhody:**
  - Silové pole je mnohem rychlejší řádné kvantové metody.
  - Reprezetnace pomocí gausovských balíčků je také a názorn, a dovoluje intuitivně pochipit and demonstrovat kvantovou chemiii
- **Nevýhody:**
  - V současné podobně nemůže přesnost konkurovat pravým kvantovým metodám v chemii. Je ale dostatečná pro popis chování látek za vysokých teplot (např. v plasmatu), nebo v rychlých procesech kdy interaguje vyskoenergetický elektron s látkout.
  - Klíčovým zjednodušením je předpoklad kulově symetrické gaussovské vlnové funkce elektronu, a popis *Pauliho repulze* mezi elektrony, která vychází z toho, že elektrony jsou fermiony a nemohou obsadit stejný kvantový stav.

#### Pauliho repulze

- V kvantovém popisu nejde o skutečnou repulzi, ale o důsledek nerozlišitelnosti fermionů. Pokud by dva elektrony byly ve stejném kvantovém stavu, jednalo by se o tentýž elektron (splynuly by).
- Protože mi nahrzuje kvantové částice kasickými, simulujeme jejich nerozlišitelnost odpudivou interakcí mezi elektrony.
- Pauliho repulze se odhaduje pomocí ortogonalizační energie mezi vlnovými funkcemi elektronů, ale tento odhad je velmi citlivý na tvar vlnových funkcí. Gaussový předpoklad vede k velké chybě.
- Podobně je tomu i u interakce valenčních elektronů s jádrem, kde se uplatňuje vliv stínění vnítřními slupkami (valenční elektron neinteraguje jen s jadrem ale taky s elektrony ve vnítřních slupkach).

## Cíl projektu

Cílem projektu je nalézt vhodnou opravu funkcí popisujících:

1. Interakci valenčních elektronů se stíněným jádrem $E_{ae}(r_{ae}, s_e)$.
2. Interakci mezi dvěma elektrony $E_{ij}(r_{ij}, s_i, s_j)$.

Tyto opravy by měly zajistit, že:

- Molekuly budou relaxovat do správné geometrie.
- Závislost celkove energie $E(\vec r)$ na poloze atomu v molekule $\vec r$ (jinak nazyvana Potenciální energetická plocha (PES)) bude co nejpřesněji kopírovat výsledky plnohodnotných kvantových metod (např. teorii funkcionálu hustoty, DFT).

### Jak nato

- Aby jsme mohli systematicky hledat nejlepší tvar funckí $E_{ae}(r_{ae}, s_e)$ a $E_{ij}(r_{ij}, s_i, s_j)$ , resp. jejich oprav, musíme je nějak vhodně reprezentovat, tak aby závisely na konečném počtu parametrů (čísel), které budeme měnit (např. náhodně), a budeme hledat pro které parametry $\vec{p}$ dostáváme nelepsí shodu (tedy nejmenší součet rozdílů $E_{EFF}(\vec r_k, \vec p) - E_{DFT}(\vec r_k)$ 0 pro všechny možné **zkušební (trénovací) geomtrie molekul $\vec r_k$**.
- To znamená že budeme hledat paramety $\vec p$ pro která má minimum výraz:

$$
min_{\vec p } \sum_k ( E_{EFF}(\vec r_k, \vec p ) - E_{DFT}(\vec r_k) )
$$

- Jak přesně budeme reprezetnovat $E_{ae}(r_{ae}, s_e)$ a $E_{ij}(r_{ij}, s_i, s_j)$ v závislosti na parametrech $\vec p$ **musíme teprve vymyslet**. Možnosti jsou např. součet gaussových funkcí, polynomů, reciprokých funckí,  spline ...
- Aby jsme vybrali vhodný matematický předpis pro reprezentaci opravy těchto funkcí potřebujeme nejdřívě vědět jak vypadají (jejich tvar, graft). 

## Problémy a výzvy

### 1. Skryté parametry ($R_i, s_i$)

- Polohy $R_i$ a poloměry $s_i$ elektronů nelze přímo porovnat s kvantovými výsledky, protože elektrony jsou v kvantových metodách delokalizované (typicky přes více atomů) a nemají tvar gaussianů s jedním jasným středem.
- Proto je nutné s $R_i, s_i$ pracovat jako se skrytými parametry, což komplikuje použití tradičních metod fitování.
- Plánujeme použít Monte Carlo, genetických algoritmů, nebo Swarm optimization pro hledání optimálních parametrů $\vec p$. 
   - Tyto metody parametry upravuj9 náhodně v okolí dosavadně nejlepšího řešení. 
   - Tím se také zbavíme problematického odvození derivací energie podle všech hledaných parametrů, což je jinak potřeba pro jiné fitovací methody.

### 2. Regularizace

- Regularizační vazby mohou penalizovat polohy elektronů které jsou v rosporu s chemickou intuicí (např. elektron mimo střed sigma-vazby vazby):
  - **Směry elektronů:** Znalosti z organické chemie (např. $sp$, $sp^2$, $sp^3$ hybridizace) lze využít k penalizaci nepravděpodobných poloh elektronů (tetraer, trojuhelník).
  - **Poloměry elektronů:** Mohou být odhadnuty z elektronové hustoty vypočítáné pomocí DFT.

### 3. Reprezentace funkcí

Je třeba zvolit vhodnou reprezentaci funkcí $E_{ae}(r_{ae}, s_e)$ a $E_{ij}(r_{ij}, s_i, s_j)$, aby byly dostatečně univerzální a přesné.

## Postup

1. **Generování dat:**
   - Vytvoříme velké množství konfigurací malých molekul ($\vec{r}_k$).
   - Pro každou konfiguraci vypočítáme pomocí DFT:
     - Energie.
     - Síly na atomy.
     - Atomové náboje (případně dipóly).
     - Elektronovou hustotu.

2. **Volba reprezentace funkcí:**
   - $E_{ae}(r_{ae}, s_e)$ – interakce elektronu s jádrem.
   - $E_{ij}(r_{ij}, s_i, s_j)$ – interakce dvou elektronů.

3. **Analýza geometrie:**
   - Zrelaxujeme $R_i$ a $s_i$ pro každou trenovací geometrii $\vec{r}_k$.
   - Vykreslíme histogramy vzorkovaných geometrií a určíme, kde je potřeba hustější vzorkování.

## Krátkodobý plán

- Zkusit nakonfigurovat / zkompilovat v Linuxu (Windows?)
- Rozmyslet si jestli programovat v C++ nebo v C#
   - C++ je rychlejší, vhodnější pro vědecké výpočty na superpočítači
     - bohužel se špatně hledají chyby (zápis/čtení mimo rozsah pole, segmentation fault, overflow)
   - C# je příjemnější na grafiku a GUI
     - Lepší chybové hlášky, automatická správa paměti(Garbage collector).
- Vykreslit $E_{ae}(r_{ae}, s_e)$ v prostoru (2D).
   - Matplotlib, OpenGL ?
- Napsat jednoduchý Monte Carlo optimalizační algoritmus v Pythonu.


## Technicke detaily

### Metrika poloh elektronů

Jak už jsem řekl, problém porovnávání molekulárních struktur z EFF a referečními strukturami z DFT jsou skryté parametry (konkrétně polohy $r_i$ a poloměry $s_i$ elektronů v DFT nejsou). Otázkou tedy je jak je porovnat.

#### Fitování elektronů

Intuitivně (zdánlivně) nejednoduží přístup je nejdříve z DFT výpočtu nějak odhadnout polohy a poloměry lokalizovaných elektronů ($r'_i, a $s'_i$ $) a pak jednoduše spočítat root-mean square error. 
$$RMSE = \sum_i ( |r_i - r'_i|^2 + (s_i - s'_i)^2 ) | )$$
Tento přístup je v principu mořný ale naráží na dvě komplikace:
* Hlavní problém smazdřejmě je to jak z DFT výpočtu polohy elektronů extrahovat. Existuje na to několik metod, žádná z nich ale není úplně jednoduchá, například:
   * To je možné udělat pomocí [Natural Bodn Orbitals](https://en.wikipedia.org/wiki/Natural_bond_orbital) analyzy, která se používá v članku [Modeling Electronic Response Properties with an Explicit-Electron Machine Learning Potential](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00978)
   * Další možnost je fitovat elektronovou hustotu pomocí gaussianů, což se v DFT výpočtech běžně dělá pro zrychlení výpočtu, viz např- [Fast periodic Gaussian density fitting by range separation](https://pubs.aip.org/aip/jcp/article/154/13/131104/1013204/Fast-periodic-Gaussian-density-fitting-by-range)
* Další problém je jak elektrony očíslovat (jak přiřadit indexy $i$)? Výsledná penalty function (chyba) by měla být nezávislá k permutaci pořadí elektronů.
   * Můžeme jednoduše spočítat minimum ze všem možných párů  
     $$R_RMSE = \sum_i min_j( |r_i - r'_j|^2| )$$
    ale to je jednak už výpočetně náročnější ( $O(n^2)$ ) 
* Problém je ale také to že tato metrika úplně nereflektuje jak se s polohou elektronů mění fyzikální vlastnosti. Například, pro elektrony s velkým poloměrem daleko od jádra nezákleží tolik kde jsou přesně, zatímco u lokalizovýnch elektronů blízko u jádra i malá změna polohy výrazně změní valastnosti. Ve skutečnosti totiž energie závisí na celkové elektronové hustotě která je součtem všech elektornů, nikoli na polohách jednotlivých elektronů.

#### Samplování Hustory

Proto by bylo mnohem výhodnější porovnávat přímo elektronvou hustotu, která jak název napovídá je z DFT snadno dostupná, a v kvantové mechanice má (na rozdíl od poloh jednotlivých elektronů) jasný výzmam (definuje energii a další vlastnosti). Nejjednoduží způsob je jednoduše vyjádřit hustotu v nějakých vybraných bodech v prostoru $R_j$. A pak jednoduše spočítat $$RMSE = \sum_i ( \rho_{EFF}(R_j) - \rho_{DFT}(R_j)' )^2$$. Referenční hustotu $\rho_{DFT}(R_i)$ je možné v kvantově mechanických programech snadno vypočíst. EFF hustotu $\rho_{EFF}(R_j)$ je možné v EFF vypočíst také snadno jako součet přes všechny jednoelektronové gaussiany
$\rho_{EFF}(R_j) = sum_i exp(- (|r_i-R_j|/s_i)^2 )$

Tato metoda je velice přímočará a robustní, jedný problém je že pro dobré porovnání může být potřeba velké množství bodů, které by měly být v prostoru nějak rozumně rovnoměrně rozmístěny.
  * Nejjednoduží by mohlo být rozmístit body jednoduše na 3D mřížku. To ale často není ideální. Pokud jsou boddy rozmísteny málo hustě, je taková metrika nepřesná (především proto že neměří efektivně hustotu okolo jádra a v chemických vazbách). Jednoduché řešení prostě jen zvýšit hustotu mřížky může být výpočetně náročné (rychle se dostatnem na tisíce a miliony vzorkovacích bodů)
  * Efektivní řešení by bylo umístit tyto samplovací do středů atomů a mezi ně. To lze udělat relativně snadno. Nicméně takové samplování může špatně popsat situaci kdy nějaký elektron "uteče pryč", tedy mimo vazby a třeba dál od molekuly (ionizovaný systém). Tyto elektrony bývají často velmi delokalizované
  * Kombinace obojího - jednoduchý a relativně přesný způsob je prostě zkombinovat oba přístupy. Můžeme mít velmi hrubou mřížku pro popis vzdálených (volných, rozprostřených) elektronů, a k tomu vzorkovací body ve středech atomů a vazeb.  

#### Multipolový rozvoj

Jasný fyzikální význam má elektrický dipol a quadrupol, které se dají snadno spočíst z DFT. Matematicky je dipol a quadrupol první a druhý [statistický moment](https://en.wikipedia.org/wiki/Moment_(mathematics)) funckce elektronové hustoty.
Dipolovy moment pro gaussovske funcke lze spočítat jako derivaci potenciálu, kde potenciál je $V(\vec{r}) = \frac{\text{erf}(\sqrt{2\alpha} r)}{r}$, tedy přesněji: 
$$V_\text{dip}(\vec{r}) \sim - \vec{p} \cdot \nabla \left( \frac{\text{erf}(\sqrt{a}r)}{r} \right)$$ 
  což dá:
  $$  V(\vec{r}) = - (\vec{p} \cdot \hat{r}) \left( \frac{2a e^{-a r^2}}{\sqrt{\pi} r} - \frac{\text{erf}(\sqrt{a}r)}{r^2} \right)$$
[viz. odvozeni v chatGPT](https://chatgpt.com/share/686cdf2d-dfa4-8003-aa96-41d68190af71)
 
## Teoreticke Pozadí

### Interakce elektronů v EFF

1. **Interakce elektron-jádro ( $E_{ae}(r_{ae}, s_e)$ ):**
   - $s_e$ je velikost elektronového oblaku (gaussovský profil).
   - $r_{ae}$ je vzdálenost mezi elektronem a jádrem.

2. **Interakce elektron-elektron ( $E_{ij}(r_{ij}, s_i, s_j)$ ):**
   - Závisí na velikostech $s_i$, $s_j$ dvou elektronů $i$ a $j$ a jejich vzdálenosti $r_{ij}$.

Obě funkce jsou redukovány o Coulombovu interakci, která je fyzikálně definována jako:

$$E^{COUL}(r_{ij}, s_i, s_j) = \text{erf} \left(\frac{r_{ij}}{s_{ij}}\right) \frac{ Q_i Q_j}{ r_{ij}}$$

kde $s_{ij}^2 = s_i^2 + s_j^2$.

Korektní interakce jsou tedy definovány jako:

- $$\Delta E_{ae} = E_{ae} - E^{COUL}_{ae},$$
- $$\Delta E_{ij} = E_{ij} - E^{COUL}_{ij}.$$

### Klíčové vlastnosti

- **Interakce elektron-jádro ($\Delta E_{ae}$):**
  - Závisí na typu atomu (nabítí jádra nebo iontu).
  - Plní roli pseudopotenciálu, který zahrnuje interakce s jádrem i vnítřními slupkami elektronů.

- **Interakce elektron-elektron ($\Delta E_{ij}$):**
  - Je univerzálnější, ale její ladění je složitější.
  - Lze případně modifikovat na základě typů blízkých atomů.

