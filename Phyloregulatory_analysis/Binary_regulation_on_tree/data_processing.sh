#GSM1579367_untreated_1.ucsc.bigWig is plus strand bw file from estaras et al in hg19

#liftover coverage to hg38:

CrossMap.py bigwig /workdir/tc643/liftOver_files/hg19ToHg38.over.chain GSM1579367_untreated_1.ucsc.bigWig GRO_h1_hg38_plus  

~/bedGraphToBigWig GRO_h1_hg38_plus.sorted.bgr /workdir/tc643/liftOver_files/hg38.chrom.sizes GRO_h1_hg38_plus.bw

CrossMap.py bigwig /workdir/tc643/liftOver_files/hg19ToHg38.over.chain GSM1579368_untreated_2.ucsc.bigWig GRO_h1_hg38_minus  


#use deeptools to rank GRO signal intensity in a -10bp to +8kb window around the 5' end of 5' and solo LTR7 
computeMatrix reference-point --referencePoint TSS -S /workdir/tc643/LTR_differential_activity/hESC/GRO_H1_Estaras2015/GRO_h1_hg38_plus.bw -R HERVH_full_bottom2 --beforeRegionStartLength 10 --afterRegionStartLength 8000 --binSize 2 --missingDataAsZero -o 3ilGROremap_noTop2.matrix.mat.gz -p max



wc -l 3ilGROremap.bed #720

head -102 3ilGROremap.bed > 3ilGROremap_1of7.bed

#1of7 is the "GRO-seq" mark shown on trees and heat maps


#KAP1 (h1) and ZNF chip (HEK):
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2067nnn/GSM2067350/suppl/GSM2067350%5FKAP1%5Fexo%5FH1%5Fpeaks%2Ebed%2Egz #kap1 h1
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2466nnn/GSM2466584/suppl/GSM2466584%5FZNF534%5Fpeaks%5Fprocessed%5Fscore%5Fsignal%5Fexo%2Ebed%2Egz #znf534 hek 
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2466nnn/GSM2466684/suppl/GSM2466684%5FZNF90%5Fpeaks%5Fprocessed%5Fscore%5Fsignal%5Fexo%2Ebed%2Egz #znf90 hek

#liftover to hg38
sed 's/^/chr/' GSM2067350_KAP1_exo_H1_peaks.bed > KAP1_exo_H1_peaks.bed
/home/tc643/liftOver KAP1_exo_H1_peaks.bed /workdir/tc643/liftOver_files/hg19ToHg38.over.chain KAP1_exo_H1_peaks_hg38.bed KAP1_exo_H1_peaks_unlifted.bed
/home/tc643/liftOver GSM2466584_ZNF534_peaks_processed_score_signal_exo.bed /workdir/tc643/liftOver_files/hg19ToHg38.over.chain ZNF534_exo_hek_peaks_hg38.bed ZNF534_exo_hek_peaks_unlifted.bed
/home/tc643/liftOver GSM2466684_ZNF90_peaks_processed_score_signal_exo.bed /workdir/tc643/liftOver_files/hg19ToHg38.over.chain ZNF90_exo_hek_peaks_hg38.bed ZNF90_exo_hek_peaks_unlifted.bed
