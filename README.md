RNA-GUI
=======


See programm on kirjutatud projekti jaoks, milles uuritakse toksiinide lõikekohti ribosomaalsel RNA-l.
Tegemist on graafilise liidesega meie RNA-seq andmete analüüsi tulemuste visualiseerimiseks.

Vajadus sellise programmi jaoks tulenes kahest põhjusest:

Esiteks on RNA-seq andmete maht väga suur(kümneid faile, millest igaüks mõni GB suur ja sisaldab miljoneid kirjeid) ning
nende jooksvalt töötlemine joonestamiseks on liiga ajakulukas, et olla praktiline.
Teiseks on meie uurimisrühmas ka inimesi, kes ei tegele igapäevaselt konsoolil programmide jooksutamisega ja kellel ei ole
Python arvutis.

Esimest probleemi proovisime me esmalt lahendada analüüsides läbi meie sekveneerimise andmed ning koostades neist .csv formaadis tabeleid, mida programm jooniste tegemiseks läbi vaatas. Selline käitlemine viis joonestamise aja umbes poole minutini, mis on parem kui pool päeva ent endiselt mitte praktiline. Praeguse lahendusena kasutame me .hdf5 formaati, millesse on array-de kujul salvestatud töödeldult need andmed, mida on joonise jaoks vaja. hdf formaadil on hea access time ja samuti on nii valitud andmete puhul võimalik failisuurus väga väikseks viia. Kui indekstabelite puhul oli meil endiselt ca 25x500MB faili, siis kõikide töötluste jaoks vajalik hdf on ligikaudu 6MB suur ja jooniseid kõikvõimalike võrdluste kohta erinevate andmehulkade vahel saab genereerida sekunditega.

Teise probleemi lahendavad nupud ja textboxid. Programm võimaldab teha mõne sekundiga jooniseid kõrvutades suvalisi töötlusi jaRNA lõikude otsi ning kuvada nende hulkasid nii protsentides kui absolutarvudes. Samuti näitab see joonisel positsioonile klikkides RNA nii absoluutset kui suhtelist lõikekohtade hulka selles punktis.

Programm on käigus olev töö ning juba praegu on kavas mitu funktsiooni, mida sellele lisada. Samuti on plaanis GUI-sse 
inkorporeeida kogu meie eelnev pipeline, et seda oleks võimalik kasutada tulevikus ka uute toorandmete puhul. Hetkel töötab
see vaid eelnevalt valmis tehtud .hdf5 failivormingus andmetega.
