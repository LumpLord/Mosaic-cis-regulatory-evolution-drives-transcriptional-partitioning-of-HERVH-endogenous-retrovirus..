#use deeptools to visualize reads in and around LTR7 insertions for repression marks and SOX2:
#see data origins and processing steps in parent folders (../../External_data.txt)
#see LTR7 bed files in partent folders. You could also generate your own bedfiles and achieve the same results

#LTR7B/C/Y are from old designations from RM library
#LTR7 subdivisions use 1:1 same insertions as in fig 1. Color matched
#insertions are ordered from most to least WGBS (rep 1) signal top to bottom

computeMatrix reference-point --referencePoint TSS -S WGBS_r1.bigWig WGBS_r2.bigWig GRO_h1_hg38_plus.bw k9me3.bigWig SOX2_primed_WIBR1_treat_pileup.bw -R OG_LTR7B.bed OG_LTR7C.bed OG_LTR7Y.bed 7up1.bed.txt rust_7up2.bed.txt peachgold_7u1.bed.txt brightgreen_7u2.bed.txt pink_7d1.bed.txt purples_7d2.bed.txt blueupper_7bc.bed.txt bluelower_7o.bed.txt --beforeRegionStartLength 500 --afterRegionStartLength 8000 --binSize 10 --missingDataAsZero -o LTR7sub_wgbs_plussox.matrix.mat.gz -p max 

plotHeatmap -m LTR7sub_wgbs_plussox.matrix.mat.gz -out LTR7sub_wgbsort_plussox.pdf  --colorMap Blues --zMin 0 --yMin 0 --sortUsingSamples 1 --perGroup


#use deeptools to visualize SOX2 reads in and around LTR7up and u1 insertions:
#insertions are ordered based on strength of sequence block 2b reverse strand match to SOX3 motif. There was overall few that had statistically significant hits to the HOMER SOX2 motif.
#Below from JDC:
# This file will contain the code to run homer2 find for the Sox2/3 motif on all of the LTR elements that Tommy provided 
# you in the LTR coords file found below. Specifically, I believe he wanted you to look into block 2b only, which requires
# the various scripts you made before to split the aligned sequences into various blocks based on the consensus sequence
# alignments. From there you can run homer2 find no problem! Good luck

# Download the list of Solo LTRs as well as 5' LTRs from Tommy. Based on the file info, 2ltr means 5' and Sltr means solo?
cp /workdir/tc643/LTR_differential_activity/imbeault_XO/HERVH_full_clustering_txn_info/manu_TF/hg38/LTR7only_final/nowwiththerightinput/adj/tree2/subgroups/new_groups/alltogethernow/RM/RMfullrun/crossmatch/LTR/LTR7allsubs.bed .

# Now manipulate the file so that it has all the info you need for bedtools getfasta. Also save the file for later usage.
### NOTE: You changed LTR7C1 and LTR7C2 to be just LTR7C since you don't have a consensus sequence for either of those 
### subfamilies separately. Make sure that's reasonable. Also confirm keeping those with length > 300 is also reasonable
egrep -v 'Un|_random|_alt' LTR7allsubs.bed | sed 's/YY/Y/g' | sed 's/YU2/u2/g' | sed 's/LTR7C[1-2]/LTR7C/g' | awk '{OFS="\t"; if ($4 ~ /2ltr/) print $1, $2, $3, "5UTR", $6, $4; if ($4 ~ /Sltr/) print $1, $2, $3, "SOLO", $6, $4} ' |  sed 's/_.*//g' | awk '{OFS="\t"; if ($3 - $2 > 300 && $5 == "C") print $1, $2 - 1, $3, $1 ":" $2 "-" $3 "|" $4 "_" $6, ".", "-"; if ($3 - $2 > 300 && $5 == "+") print $1, $2 - 1, $3, $1 ":" $2 "-" $3 "|" $4 "_" $6, ".", $5}' > subFams_clean.bed

# Turn the above file into .fa files corresponding to each type of subfamily via grep! 
for f in $(cut -f4 subFams_clean.bed | sed 's/.*_//g' | sort | uniq); do grep $f subFams_clean.bed | bedtools getfasta -s -name -fi /workdir/Genomes/homoSap/hg38/hg38.fa -bed stdin | sed 's/::.*)//g' > "${f}".fa; done

# Download the .fa consensus sequences you generated before from Tommy. You manually added 'LTR' to match the above file as needed
cp /workdir/jdc397/0_orthTEs/1_ChIP/6_LTR7_tommy/0_allLTR7/*_cons.fa .

# Now download the script you made before to combine the consensus sequences, align them and then split them into blocks!
cp /workdir/jdc397/0_orthTEs/1_ChIP/6_LTR7_tommy/0_allLTR7/2_tom7blocks/0_scripts.v2/LTR7_muscleConsStrand.sh LTR7_muscleConsStrandFind.sh
cp /workdir/jdc397/0_orthTEs/1_ChIP/6_LTR7_tommy/0_allLTR7/2_tom7blocks/0_scripts.v2/LTR7_blockFilt.v2.r LTR7_blockFilt.r

### NOTE: You manually changed the above scripts to handle the new jobs/names and also took the Sox3 motif from 
### /workdir/jdc397/homer/data/knownTFs/vertebrates/known.motifs called sox3Motif.motifs for use in homer2 find


#LTR7_muscleConsStrand.sh and LTR7_blockFilt.r are in this folder #TAC


#below from TAC
#generate files
#HOMERfind.beds from JDC (see accompanying code for finding and scoring SOX3 motif in block 2b)
grep "7up1" block2b_sox3_HOMERfind.bed.txt > 7up1_block2b_sox3_HOMERfind.bed
grep "7up2" block2b_sox3_HOMERfind.bed.txt > 7up2_block2b_sox3_HOMERfind.bed
grep "7u1" block2b_sox3_HOMERfind.bed.txt > 7u1_block2b_sox3_HOMERfind.bed

intersectBed -u -a <(awk '{print $1 "\t" $2 "\t" $3}' /workdir/tc643/LTR_differential_activity/imbeault_XO/HERVH_full_clustering_txn_info/manu_TF/hg38/LTR7only_final/nowwiththerightinput/adj/tree2/subgroups/new_groups/7up1.bed.txt) -b 7up1_block2b_sox3_HOMERfind.bed > 7up1_sox.bed
intersectBed -v -a <(awk '{print $1 "\t" $2 "\t" $3}' /workdir/tc643/LTR_differential_activity/imbeault_XO/HERVH_full_clustering_txn_info/manu_TF/hg38/LTR7only_final/nowwiththerightinput/adj/tree2/subgroups/new_groups/7up1.bed.txt) -b 7up1_block2b_sox3_HOMERfind.bed > 7up1_nosox.bed

intersectBed -u -a <(awk '{print $1 "\t" $2 "\t" $3}' /workdir/tc643/LTR_differential_activity/imbeault_XO/HERVH_full_clustering_txn_info/manu_TF/hg38/LTR7only_final/nowwiththerightinput/adj/tree2/subgroups/new_groups/rust_7up2.bed.txt) -b 7up2_block2b_sox3_HOMERfind.bed > 7up2_sox.bed
intersectBed -v -a <(awk '{print $1 "\t" $2 "\t" $3}' /workdir/tc643/LTR_differential_activity/imbeault_XO/HERVH_full_clustering_txn_info/manu_TF/hg38/LTR7only_final/nowwiththerightinput/adj/tree2/subgroups/new_groups/rust_7up2.bed.txt) -b 7up2_block2b_sox3_HOMERfind.bed > 7up2_nosox.bed

intersectBed -u -a <(awk '{print $1 "\t" $2 "\t" $3}' /workdir/tc643/LTR_differential_activity/imbeault_XO/HERVH_full_clustering_txn_info/manu_TF/hg38/LTR7only_final/nowwiththerightinput/adj/tree2/subgroups/new_groups/peachgold_7u1.bed.txt) -b 7u1_block2b_sox3_HOMERfind.bed > 7u1_sox.bed
intersectBed -v -a <(awk '{print $1 "\t" $2 "\t" $3}' /workdir/tc643/LTR_differential_activity/imbeault_XO/HERVH_full_clustering_txn_info/manu_TF/hg38/LTR7only_final/nowwiththerightinput/adj/tree2/subgroups/new_groups/peachgold_7u1.bed.txt) -b 7u1_block2b_sox3_HOMERfind.bed > 7u1_nosox.bed

cat 7up1_sox.bed 7up1_nosox.bed > 7up1_soxsort.bed #insertions with sox scores are ontop in correct otientation
cat 7up2_sox.bed 7up2_nosox.bed > 7up2_soxsort.bed
cat 7u1_sox.bed 7u1_nosox.bed > 7u1_soxsort.bed


computeMatrix reference-point --referencePoint TSS -S /workdir/tc643/LTR_differential_activity/hESC/GRO_H1_Estaras2015/GRO_h1_hg38_plus.bw SOX2_primed_WIBR1_treat_pileup.bw -R 7up1_soxsort.bed 7up2_soxsort.bed 7u1_soxsort.bed --beforeRegionStartLength 500 --afterRegionStartLength 8000 --binSize 10 --missingDataAsZero -o LTR7up_7u_sox3motif.matrix.mat.gz -p max 


plotHeatmap -m LTR7up_7u_sox3motif.matrix.mat.gz -out LTR7up_7u_sox3motif.pdf  --colorMap Blues --zMin 0 --yMin 0 --sortRegions no --perGroup

#######
#to test if there was more bw signal on LTR7 with SOX3 motif in lbock 2a vs those with none:
#######

multiBigwigSummary BED-file \
 --bwfiles SOX2_primed_WIBR1_treat_pileup.bw \
 --BED 7up1_sox.bed \
 --labels sox2 \
 -out 7up1_sox3motif_sox2signal.npz --outRawCounts 7up1_sox3motif_sox2signal.tab

multiBigwigSummary BED-file \
 --bwfiles SOX2_primed_WIBR1_treat_pileup.bw \
 --BED 7up1_nosox.bed \
 --labels sox2 \
 -out 7up1_nosox3motif_sox2signal.npz --outRawCounts 7up1_nosox3motif_sox2signal.tab

awk '{print $4}' 7up1_nosox3motif_sox2signal.tab #avg=4.19
awk '{print $4}' 7up1_sox3motif_sox2signal.tab #avg=6.23
#mann whitney p=.04136

multiBigwigSummary BED-file \
 --bwfiles SOX2_primed_WIBR1_treat_pileup.bw \
 --BED 7up2_sox.bed \
 --labels sox2 \
 -out 7up2_sox3motif_sox2signal.npz --outRawCounts 7up2_sox3motif_sox2signal.tab

multiBigwigSummary BED-file \
 --bwfiles SOX2_primed_WIBR1_treat_pileup.bw \
 --BED 7up2_nosox.bed \
 --labels sox2 \
 -out 7up2_nosox3motif_sox2signal.npz --outRawCounts 7up2_nosox3motif_sox2signal.tab

awk '{print $4}' 7up2_nosox3motif_sox2signal.tab #avg=5.55
awk '{print $4}' 7up2_sox3motif_sox2signal.tab #avg=8.67
# mann whitney p=0.05744

multiBigwigSummary BED-file \
 --bwfiles SOX2_primed_WIBR1_treat_pileup.bw \
 --BED 7u1_sox.bed \
 --labels sox2 \
 -out 7u1_sox3motif_sox2signal.npz --outRawCounts 7u1_sox3motif_sox2signal.tab

multiBigwigSummary BED-file \
 --bwfiles SOX2_primed_WIBR1_treat_pileup.bw \
 --BED 7u1_nosox.bed \
 --labels sox2 \
 -out 7u1_nosox3motif_sox2signal.npz --outRawCounts 7u1_nosox3motif_sox2signal.tab


awk '{print $4}' 7u1_nosox3motif_sox2signal.tab #avg=5.55
awk '{print $4}' 7u1_sox3motif_sox2signal.tab #avg=8.67

#p=0.04136

#no p survives MTC 