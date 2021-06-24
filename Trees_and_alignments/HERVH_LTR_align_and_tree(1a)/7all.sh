#7all_350.fa is all 5' and solo HERVH LTR from remasking >350bp (see methods, files folder)


/programs/mafft/bin/mafft --thread 100 --auto 7all_350.fa > 7all_350_mafft.fa

Strategy:
 FFT-NS-2 (Fast but rough)
 Progressive method (guide trees were built 2 times.)

sed -i 's/:/_/g' 7all_350_mafft.fa
sed -i 's/)/_/g' 7all_350_mafft.fa
sed -i 's/(/_/g' 7all_350_mafft.fa
sed -i 's/;/_/g' 7all_350_mafft.fa
sed -i 's/,/_/g' 7all_350_mafft.fa



export LD_LIBRARY_PATH=/usr/local/gcc-7.3.0/lib:/usr/local/gcc-7.3.0/lib64
/programs/prank-170427/bin/prank -d=7all_350_mafft.fa -o=./7all_350_mafft_prank_F1 -showanc -support -njtree -uselogs -prunetree -prunedata -F -showevents


#trimal
/home/tc643/trimal/source/trimal -in 7all_350_mafft_prank_F1.best.fas -out 7all_350_mafft_prank_F1.best.trimal.fa -gt .01 -fasta




/programs/iqtree-1.6.10-Linux/bin/iqtree  -s  7all_350_mafft_prank_F1.best.trimal.fa -nt AUTO -pre 7all_350_mafft_prank_F1.best.trimal_iqtree_input -m MFP -bb 6000 -asr -ntmax 100 -minsup .95


