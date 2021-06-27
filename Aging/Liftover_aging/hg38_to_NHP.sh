#LTR7_custom_subs.bed is a simple .bed with group definitions from figure 1B (LTR7) 
#USCS liftover chains have original file names
#pt6 = panTro6
#gg5 = gorGor5
#pa3 = ponAbe3
#NL3 = nomLeu3
#rm10 = rheMac10

#combine new LTR7 solo and 5' annotations (fig 1b) with previously annotated LTR7B/C/Y (RepeatMasker version 4.0.5 repeat Library 20140131 - "hg38.out" here)



#get LTR7B/C/Y RM coordinates:
grep -w 'LTR7B' /workdir/tc643/RMout/hg38.out | awk '{print $5 "\t" $6 "\t" $7 "\t" $9 "\t" "LTR7B"}' > LTR7B.bed
grep -w 'LTR7C' /workdir/tc643/RMout/hg38.out | awk '{print $5 "\t" $6 "\t" $7 "\t" $9 "\t" "LTR7C"}' > LTR7C.bed
grep -w 'LTR7Y' /workdir/tc643/RMout/hg38.out | awk '{print $5 "\t" $6 "\t" $7 "\t" $9 "\t" "LTR7Y"}' > LTR7Y.bed

#using above to extract only LTR from OneCode run (see "HERVH_OneCode_dfam_runs.sh")
intersectBed -u -wa -a /workdir/tc643/RMout/hg38_oneCode/7B/LTR7B_full.bed -b LTR7B.bed > LTR7B_full_LTRs.bed
grep "+" LTR7B_full_LTRs.bed | awk '{print $1 "\t" $2 "\t" $3-2000 "\t" $4 "\t" $5 "\t" $6}' > LTR7B_full_LTRs_plushalf.bed
grep -v "+" LTR7B_full_LTRs.bed | awk '{print $1 "\t" $2+2000 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > LTR7B_full_LTRs_minushalf.bed
cat LTR7B_full_LTRs_*half.bed > LTR7B_full_LTRshalves.bed
#chr9    11103109        11102669 and chr12   54933682        54932632 causing issue due to being very short (not actually full length) - removed. 
intersectBed -a <(sort -k1,1 -k2,2n LTR7B_full_LTRshalves.bed) -b <(sort -k1,1 -k2,2n LTR7B.bed) > LTR7B_full_LTRs.bed

intersectBed -u -wa -a /workdir/tc643/RMout/hg38_oneCode/7B/LTR7B_solo.bed -b LTR7B.bed > LTR7B_solo_LTRs.bed
cat LTR7B_full_LTRs.bed LTR7B_solo_LTRs.bed > LTR7B_ALL.bed





intersectBed -u -wa -a /workdir/tc643/RMout/hg38_oneCode/7C/LTR7C_full.bed -b LTR7C.bed > LTR7C_full_LTRs.bed
grep "+" LTR7C_full_LTRs.bed | awk '{print $1 "\t" $2 "\t" $3-2000 "\t" $4 "\t" $5 "\t" $6}' > LTR7C_full_LTRs_plushalf.bed
grep -v "+" LTR7C_full_LTRs.bed | awk '{print $1 "\t" $2+2000 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > LTR7C_full_LTRs_minushalf.bed
cat LTR7C_full_LTRs_*half.bed > LTR7C_full_LTRshalves.bed
#chr4    124878070       124877475 causing issue due to being very short (not actually full length) - removed. 
intersectBed -a <(sort -k1,1 -k2,2n LTR7C_full_LTRshalves.bed) -b <(sort -k1,1 -k2,2n LTR7C.bed) > LTR7C_full_LTRs.bed

intersectBed -u -wa -a /workdir/tc643/RMout/hg38_oneCode/7C/LTR7C_solo.bed -b LTR7C.bed > LTR7C_solo_LTRs.bed

cat LTR7C_full_LTRs.bed LTR7C_solo_LTRs.bed > LTR7C_ALL.bed




intersectBed -u -wa -a /workdir/tc643/RMout/hg38_oneCode/7C/LTR7Y_full.bed -b LTR7Y.bed > LTR7Y_full_LTRs.bed
grep "+" LTR7Y_full_LTRs.bed | awk '{print $1 "\t" $2 "\t" $3-2000 "\t" $4 "\t" $5 "\t" $6}' > LTR7Y_full_LTRs_plushalf.bed
grep -v "+" LTR7Y_full_LTRs.bed | awk '{print $1 "\t" $2+2000 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > LTR7Y_full_LTRs_minushalf.bed
cat LTR7Y_full_LTRs_*half.bed > LTR7Y_full_LTRshalves.bed
#chr15   55717804        55716624 chr2    57192263        57191496 chr3    156047103       156045445 causing issue due to being very short (not actually full length) - removed. 
intersectBed -a <(sort -k1,1 -k2,2n LTR7Y_full_LTRshalves.bed) -b <(sort -k1,1 -k2,2n LTR7Y.bed) > LTR7Y_full_LTRs.bed


intersectBed -u -wa -a /workdir/tc643/RMout/hg38_oneCode/7C/LTR7Y_solo.bed -b LTR7Y.bed > LTR7Y_solo_LTRs.bed

cat LTR7Y_full_LTRs.bed LTR7Y_solo_LTRs.bed > LTR7Y_ALL.bed



#combine LTR7B/C/Y LTR coordinates with new definitions for LTR7 (fig1b)
cat LTR7_custom_subs.bed /workdir/tc643/LTR_differential_activity/other_subfamilies/LTR7C_ALL.bed  /workdir/tc643/LTR_differential_activity/other_subfamilies/LTR7B_ALL.bed  /workdir/tc643/LTR_differential_activity/other_subfamilies/LTR7Y_ALL.bed > all_named.bed


#######################
#LTR7 subgroup liftover:
#LiftOver human to chimp
liftOver ../../all_named.bed /workdir/tc643/liftOver_files/hg38ToPanTro6.over.chain LTR7subs_all_to_pt6.bed LTR7subs_all_notto_pt6.bed
#ReciprocalLO chimp to human
liftOver LTR7subs_all_to_pt6.bed /workdir/tc643/liftOver_files/panTro6ToHg38.over.chain LTR7subs_all_RLO_pt6.bed LTR7subs_all_noRLO_pt6.bed

#LO human to GOrilla
liftOver ../../all_named.bed /workdir/tc643/liftOver_files/hg38ToGorGor5.over.chain LTR7subs_all_to_gg5.bed LTR7subs_all_notto_gg5.bed
#RLO GO to human
liftOver LTR7subs_all_to_gg5.bed /workdir/tc643/liftOver_files/gorGor5ToHg38.over.chain LTR7subs_all_RLO_gg5.bed LTR7subs_all_noRLO_gg5.bed

#LO human to ORangutan
liftOver ../../all_named.bed /workdir/tc643/liftOver_files/hg38ToPonAbe3.over.chain LTR7subs_all_to_pa3.bed LTR7subs_all_notto_pa3.bed
#RLO OR to human
liftOver LTR7subs_all_to_pa3.bed /workdir/tc643/liftOver_files/ponAbe3ToHg38.over.chain LTR7subs_all_RLO_pa3.bed LTR7subs_all_noRLO_pa3.bed

#LO human to gibbon
liftOver ../../all_named.bed /workdir/tc643/liftOver_files/hg38ToNomLeu3.over.chain LTR7subs_all_to_NL3.bed LTR7subs_all_notto_NL3.bed
#RLO gi to human
liftOver LTR7subs_all_to_NL3.bed /workdir/tc643/liftOver_files/nomLeu3ToHg38.over.chain LTR7subs_all_RLO_NL3.bed LTR7subs_all_noRLO_NL3.bed

#LO human to RhesusMacaque
liftOver ../../all_named.bed /workdir/tc643/liftOver_files/hg38ToRheMac10.over.chain LTR7subs_all_to_RM10.bed LTR7subs_all_notto_rm10.bed
#RLO RM to human
liftOver LTR7subs_all_to_RM10.bed /workdir/tc643/liftOver_files/rheMac10ToHg38.over.chain LTR7subs_all_RLO_RM10.bed LTR7subs_all_noRLO_rm10.bed


