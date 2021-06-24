#LTR7only_sort.fa is 5' and solo LTR7 >350bp (see methods, files folder)

/programs/mafft/bin/mafft --thread 100 --auto LTR7only_sort.fa > LTR7only_clean_mafft.fa
Strategy:
 FFT-NS-2 (Fast but rough)
 Progressive method (guide trees were built 2 times.)

sed -i 's/:/_/g' LTR7only_clean_mafft.fa
sed -i 's/)/_/g' LTR7only_clean_mafft.fa
sed -i 's/(/_/g' LTR7only_clean_mafft.fa
sed -i 's/;/_/g' LTR7only_clean_mafft.fa
sed -i 's/,/_/g' LTR7only_clean_mafft.fa

/programs/prank-170427/bin/prank -d=../LTR7only_clean_mafft.fa -o=./LTR7only_clean_mafft_prank_F1 -showanc -support -njtree -uselogs -prunetree -prunedata -F
#pruned outliers: 
#2 34800999 34806725
#y 14465541 14465939
#y 15098843 15103998
#y 21044747 21045245


/programs/iqtree-1.6.10-Linux/bin/iqtree -s  LTR7only_clean_mafft_prank_F1.trimal_pruned.fa -nt AUTO -pre LTR7only_clean_mafft_prank_F1_pruned.trimal_iqtree_input -m MFP -bb 6000 -asr -ntmax 100 -minsup .95
#model:
Akaike Information Criterion:           TIM3+F+R8
Corrected Akaike Information Criterion: JC+G4
Bayesian Information Criterion:         TVMe+R8
Best-fit model: TVMe+R8 chosen according to BIC

