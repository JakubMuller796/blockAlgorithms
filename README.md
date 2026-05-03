# blockAlgorithms

V tomto repozitáři jsou MATLAB skripty a funkce použité k numerickým experimentům v mé bakalářské práci.

Byly implementovány a testovány čtyři blokové algoritmy:
- Blokový Lanczosův algoritmus (Golub, Underwood, 1977) 
- BCG (Dianne O'Leary, 1980),  
- DRBCG (Tichý, Meraunt, Šimonová, 2025),  
- BFBCG (Yi, Li, 2017),  

tak, aby bylo možné mé výsledky zreprodukovat.

## Nastavení

Doporučuje se nastavit pracovní adresář na '.../blockAlgorithms', protože výsledky se standardně ukládají do složky '/exp' (pro velké matice mohou být výpočty dlouhé, proto je vhodné výsledky ukládat namísto neustálého přepočítávání).

Samotné algoritmy, vykreslovací funkce a pomocné funkce se nacházejí ve složce '/src'.

Ve složce '/exp' jsou předpočítané výsledky, tedy je možné spustit pouze vykreslovací části skriptů (po načtení proměnné 'matrix').

Složka '/matrix' obsahuje matice použité v experimentech a dvě menší testovací matice.

Kódy jsou napsány tak, aby je bylo možné použít s libovolnými maticemi uloženými do složky '/matrix' a specifikovanými v příslušném skriptu. Výchozí nastavení v jednotlivých skriptech odpovídá experimentům z mé bakalářské práce, lze je nicméně snadno upravit.

---

## Parametry

- 'm' – počet pravých stran  
- 'tol' – tolerance pro konvergenci
- 'rankrevtol' – tolerance pro odhad numerické hodnosti (pro '0' se numerická hodnost nezjišťuje)  
- 'maxit' – maximální počet iterací  
- 'coef' (0/1) – výpočet tridiagonální matice  
- 'reo' (0/1/...) – vícenásobná reortogonalizace  
- 'tillConv' (0/1) – pro '1' se algoritmy zastaví po dosažení tolerance 'tol'; pro '0' provedou plný počet iterací 'maxit'  

---

## Návratové hodnoty

Algoritmy BCG, DRBCG a BFBCG vracejí hodnoty:

- 'x' — vypočtená aproximace řešení  
- 'T' — bloková tridiagonální matice z blokového Lanczosova algoritmu (pouze BCG a DRBCG)  
- 'omega' — analogie relativní A-normy počítané přímo v algoritmech; pouze pokud je zadáno 'xex' přesné řešení
- 'numrank' — pro 'rankrevtol > 0' je numerická hodnost odhadnuta v každé iteraci pomocí SVD, v případě změny se uloží numerická hodnost a iterace (výpočty pak ale mohou trvat déle)  
- 'conv' — počet iterací a čas do zkonvergování (vzhledem k zadané toleranci 'tol')

Blokový Lanczosův algoritmus vrací hodnoty:

- 'V' — matice s (teoreticky) ortonormálními blokovými vektory  
- 'T' — bloková tridiagonální matice  
- 'numrank' — stejné jako výše  

---

## Hlavní skripty

- `exp1_fitting.m`  
  První experiment: porovnání blokové tridiagonální matice počítané blokovým Lanczosovým algoritmem, BCG a DRBCG.  
  Ztratí-li blokové vektory v blokovém Lanczosově algoritmu hodnost, pak nejsou zaručeny správné výsledky.

- `exp2_compar.m`  
  Druhý experiment: porovnání algoritmů BCG, DRBCG a BFBCG pro čtyři různé hodnoty 'm'.

- `exp3_defl.m`  
  Třetí experiment: sledování numerické sloupcové hodnosti v závislosti na různých hodnotách 'rankrevtol'.  
  Funguje i pro více tolerancí, stačí je přidat do vektoru 'rankrevtols' a doplnit odpovídající počet hodnot 'svd' nebo 'qr' do vektoru 'facs'.

- `xxCG_testing.m`  
  Testovací skript určený pro experimentování pouze s jedním algoritmem.

---

## Poznámky

Pro teoretické a experimentální účely jsou ukládány některé další hodnoty a je možné upravovat více parametrů (například rozklady v algoritmech DRBCG a BFBCG).

Ve druhých sekcích lze také nastavit vykreslování více informací v grafech.

Vše by mělo být plně funkční (případně snadno opravitelné) a připravené pro různé experimenty.
