#repeatmasking hg38 (alt removed) with custom library
#

/home/tc643/RepeatMasker-4.1.1/RepeatMasker/RepeatMasker -e crossmatch -lib LTR7cons_bcydfam_int_mismatches.fa -pa 16 -a -s -no_is -dir /workdir/tc643/LTR_differential_activity/imbeault_XO/HERVH_full_clustering_txn_info/manu_TF/hg38/LTR7only_final/nowwiththerightinput/adj/tree2/subgroups/new_groups/alltogethernow/RM/RMfullrun/crossmatch /workdir/Genomes/homoSap/hg38/hg38.fa

#rm alt chromosomes
grep -v "alt" hg38.fa.out > hg38.fa_notalt.out
