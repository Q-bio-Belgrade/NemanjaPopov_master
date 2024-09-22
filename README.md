# MATLAB kod za reprodukovanje rezultata

Ovaj repozitorijum sadrži MATLAB kod za reprodukovanje rezultata iz master rada:

**"Inverzna optimizacija modela CRISPR/Cas sistema tipa I-E u Escherichia coli iz simuliranih podataka"** 
od Nemanje Popov

## Skripte

Skripte u ovom repozitorijumu su obezbeđene radi reprodukcije rezultata predstavljenih u master radu.

## crispr.m

### Opis:
Ova skripta simulira CRISPR-Cas sistema za originalne vrednosti parametara i daje grafik dinamike crRNK naspram simuliranih podataka sa slike 2 iz rada.

## gradient_descent.m

### Opis:
U ovoj skripti se nalazi kod sa kojim su dobijeni rezultati algoritma opadajućeg gradijenta, slike 3-5 iz rada.

## genetic_algorithm.m

### Opis:
Pokretanjem skripte se dobijaju rezultati genetskog algoritma, slike 6-8 iz rada.

## grid.m

### Opis:
MATLAB kod za inferenciju vrednosti parametara mrežnom pretragom primenjenom u master radu.

## Beleške

- Kod je pisan u MATLAB verziji R2024a.
- Proverite da su sve MATLAB alatne kutije instalirane kako bi skripte ispravno radile. (Statistics and Machine Learning Toolbox, Optimization i Global Optimization Toolbox)
- Fajl perturbed_data_high.mat sadrži simulirane podatke na osnovu koji se radi inferencija parametara. Neophodno je držati ovaj fajl u istom folderu kao i skripte kako bi kod radio ispravno.

## Licenca

Ovaj projekat je licenciran pod MIT Licencom. Pogledajte [MIT License](https://opensource.org/licenses/MIT) za detalje.




