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

---

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

