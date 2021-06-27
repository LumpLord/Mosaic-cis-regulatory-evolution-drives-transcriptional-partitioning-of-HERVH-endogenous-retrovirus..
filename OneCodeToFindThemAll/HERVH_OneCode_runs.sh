#OneCodeToFindThemAll runs of REMASKED hg38 with LTR7 subgroup definitions.
#rename_mergedLTRelements_v6.0.pl will cat OneCode output and give Full, Solo, and Truncated counts (insertion counts, not LTR counts which you would get from RM.out)
cd ./7b
vi LTR7b.dico
#HERVH-int    LTR7B1
grep "HERVH-int\|LTR7B1" ../../hg38.fa_notalt.out > LTR7B.out

#rm "*" character and put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7B.out 
sed -i 's/\*//g' LTR7B.out 


~/one_code_to_find_them_all.pl -rm LTR7B.out -ltr LTR7B.dico > LTR7B.OC



cat *ltr.csv > LTR7b.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7b.csv -ltr LTR7B1=414 -int HERVH-int -ilen 7713 -rn LTR7B  
grep '2ltr_int' LTR7b.csv.newnames.cordinates.txt > LTR7B_full.bed
grep '1ltr_int' LTR7b.csv.newnames.cordinates.txt > LTR7B_truncated.bed
grep 'Sltr' LTR7b.csv.newnames.cordinates.txt > LTR7B_solo.bed
wc -l *7B_*.bed 
   81 LTR7B_full.bed
  367 LTR7B_solo.bed
   99 LTR7B_truncated.bed
  547 total
#699 via grep "LTR7B1" ../hg38.fa.out

#++++++++++++++++++++++++++++++++++++++++++++++

cd ../7bc
vi LTR7bc.dico
#HERVH-int    LTR7bc
grep "HERVH-int\|LTR7bc" ../../hg38.fa_notalt.out > LTR7bc.out

#rm "*" character and put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7bc.out 
#sed -i 's/\*//g' LTR7B.out 


~/one_code_to_find_them_all.pl -rm LTR7bc.out -ltr LTR7bc.dico > LTR7bc.OC



cat *ltr.csv > LTR7bc.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7bc.csv -ltr LTR7bc=395 -int HERVH-int -ilen 7713 -rn LTR7bc 
grep '2ltr_int' LTR7bc.csv.newnames.cordinates.txt > LTR7bc_full.bed
grep '1ltr_int' LTR7bc.csv.newnames.cordinates.txt > LTR7bc_truncated.bed
grep 'Sltr' LTR7bc.csv.newnames.cordinates.txt > LTR7bc_solo.bed
wc -l *7bc_*.bed 
   27 LTR7bc_full.bed
  161 LTR7bc_solo.bed
  119 LTR7bc_truncated.bed
  307 total
#340 via grep -c "LTR7bc" ../hg38.fa.out

#++++++++++++++++++++++++++++++++++++++++++++++

cd ../7c1
vi LTR7c1.dico
#HERVH-int    LTR7C1_cons
grep "HERVH-int\|LTR7C1" ../../hg38.fa_notalt.out > LTR7c1.out

#rm "*" character and put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7c1.out 
#sed -i 's/\*//g' LTR7B.out 


~/one_code_to_find_them_all.pl -rm LTR7c1.out -ltr LTR7c1.dico > LTR7c1.OC



cat *ltr.csv > LTR7c1.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7c1.csv -ltr LTR7C1_cons=440 -int HERVH-int -ilen 7713 -rn LTR7C1 
grep '2ltr_int' LTR7c1.csv.newnames.cordinates.txt > LTR7c1_full.bed
grep '1ltr_int' LTR7c1.csv.newnames.cordinates.txt > LTR7c1_truncated.bed
grep 'Sltr' LTR7c1.csv.newnames.cordinates.txt > LTR7c1_solo.bed
wc -l *7c1_*.bed 
  19 LTR7c1_full.bed
  165 LTR7c1_solo.bed
   61 LTR7c1_truncated.bed
  245 total
# via grep -c "LTR7c1" ../hg38.fa.out

#++++++++++++++++++++++++++++++++++++++++++++++

cd ../7c2
vi LTR7C2.dico
#HERVH-int    LTR7C2_cons
grep "HERVH-int\|LTR7C2" ../../hg38.fa_notalt.out > LTR7C2.out

#rm "*" character and put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7C2.out 
#sed -i 's/\*//g' LTR7B.out 


~/one_code_to_find_them_all.pl -rm LTR7C2.out -ltr LTR7C2.dico > LTR7C2.OC



cat *ltr.csv > LTR7C2.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7C2.csv -ltr LTR7C2_cons=436 -int HERVH-int -ilen 7713 -rn LTR7C2 
grep '2ltr_int' LTR7C2.csv.newnames.cordinates.txt > LTR7C2_full.bed
grep '1ltr_int' LTR7C2.csv.newnames.cordinates.txt > LTR7C2_truncated.bed
grep 'Sltr' LTR7C2.csv.newnames.cordinates.txt > LTR7C2_solo.bed
wc -l *7C2_*.bed 
  15 LTR7C2_full.bed
 123 LTR7C2_solo.bed
  29 LTR7C2_truncated.bed
 167 total
# via grep -c "LTR7C2" ../hg38.fa.out

#++++++++++++++++++++++++++++++++++++++++++++++

cd ../7d1
vi LTR7d1.dico
#HERVH-int    LTR7d1
grep "HERVH-int\|LTR7d1" ../../hg38.fa_notalt.out > LTR7d1.out

# put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7d1.out 



~/one_code_to_find_them_all.pl -rm LTR7d1.out -ltr LTR7d1.dico > LTR7d1.OC



cat *ltr.csv > LTR7d1.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7d1.csv -ltr LTR7d1=416 -int HERVH-int -ilen 7713 -rn LTR7d1 
grep '2ltr_int' LTR7d1.csv.newnames.cordinates.txt > LTR7d1_full.bed
grep '1ltr_int' LTR7d1.csv.newnames.cordinates.txt > LTR7d1_truncated.bed
grep 'Sltr' LTR7d1.csv.newnames.cordinates.txt > LTR7d1_solo.bed
wc -l *7d1_*.bed
   72 LTR7d1_full.bed
   77 LTR7d1_solo.bed
   66 LTR7d1_truncated.bed
  215 total
#310 via grep -c "LTR7d1" ../hg38.fa.out

#++++++++++++++++++++++++++++++++++++++++++++++

cd ../7d2
vi LTR7d2.dico
#HERVH-int    LTR7d2
grep "HERVH-int\|LTR7d2" ../../hg38.fa_notalt.out > LTR7d2.out

# put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7d2.out 



~/one_code_to_find_them_all.pl -rm LTR7d2.out -ltr LTR7d2.dico > LTR7d2.OC



cat *ltr.csv > LTR7d2.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7d2.csv -ltr LTR7d2=416 -int HERVH-int -ilen 7713 -rn LTR7d2 
grep '2ltr_int' LTR7d2.csv.newnames.cordinates.txt > LTR7d2_full.bed
grep '1ltr_int' LTR7d2.csv.newnames.cordinates.txt > LTR7d2_truncated.bed
grep 'Sltr' LTR7d2.csv.newnames.cordinates.txt > LTR7d2_solo.bed
wc -l *7d2_*.bed
  177 LTR7d2_full.bed
  151 LTR7d2_solo.bed
  117 LTR7d2_truncated.bed
  445 total
#659 via grep -c "LTR7d2" ../hg38.fa.out

#++++++++++++++++++++++++++++++++++++++++++++++

cd ../7u1
vi LTR7u1.dico
#HERVH-int    LTR7u1
grep "HERVH-int\|LTR7u1" ../../hg38.fa_notalt.out > LTR7u1.out

# put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7u1.out 



~/one_code_to_find_them_all.pl -rm LTR7u1.out -ltr LTR7u1.dico > LTR7u1.OC



cat *ltr.csv > LTR7u1.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7u1.csv -ltr LTR7u1=455 -int HERVH-int -ilen 7713 -rn LTR7u1 
grep '2ltr_int' LTR7u1.csv.newnames.cordinates.txt > LTR7u1_full.bed
grep '1ltr_int' LTR7u1.csv.newnames.cordinates.txt > LTR7u1_truncated.bed
grep 'Sltr' LTR7u1.csv.newnames.cordinates.txt > LTR7u1_solo.bed
wc -l *7u1_*.bed
  21 LTR7u1_full.bed
  90 LTR7u1_solo.bed
  32 LTR7u1_truncated.bed
 143 total
#174 via grep -c "LTR7u1" ../hg38.fa.out

#++++++++++++++++++++++++++++++++++++++++++++++

cd ../7u2
vi LTR7u2.dico
#HERVH-int    LTR7U2_cons
grep "HERVH-int\|LTR7YU2_cons" ../../hg38.fa_notalt.out > LTR7u2.out

# put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7u2.out 



~/one_code_to_find_them_all.pl -rm LTR7u2.out -ltr LTR7u2.dico > LTR7u2.OC



cat *ltr.csv > LTR7u2.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7u2.csv -ltr LTR7YU2_cons=433 -int HERVH-int -ilen 7713 -rn LTR7YU2_cons
grep '2ltr_int' LTR7u2.csv.newnames.cordinates.txt > LTR7u2_full.bed
grep '1ltr_int' LTR7u2.csv.newnames.cordinates.txt > LTR7u2_truncated.bed
grep 'Sltr' LTR7u2.csv.newnames.cordinates.txt > LTR7u2_solo.bed
wc -l *7u2_*.bed
   62 LTR7u2_full.bed
   78 LTR7u2_solo.bed
   55 LTR7u2_truncated.bed
  195 total
#257 via grep -c "LTR7u2" ../hg38.fa.out

#++++++++++++++++++++++++++++++++++++++++++++++

cd ../7o
vi LTR7o.dico
#HERVH-int    LTR7o_cons
grep "HERVH-int\|LTR7o" ../../hg38.fa_notalt.out > LTR7o.out

# put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7o.out 



~/one_code_to_find_them_all.pl -rm LTR7o.out -ltr LTR7o.dico > LTR7o.OC



cat *ltr.csv > LTR7o.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7o.csv -ltr LTR7o=371 -int HERVH-int -ilen 7713 -rn LTR7o
grep '2ltr_int' LTR7o.csv.newnames.cordinates.txt > LTR7o_full.bed
grep '1ltr_int' LTR7o.csv.newnames.cordinates.txt > LTR7o_truncated.bed
grep 'Sltr' LTR7o.csv.newnames.cordinates.txt > LTR7o_solo.bed
wc -l *7o_*.bed
   39 LTR7o_full.bed
  141 LTR7o_solo.bed
   63 LTR7o_truncated.bed
  243 total
#290 via grep -c "LTR7o" ../hg38.fa.out

#++++++++++++++++++++++++++++++++++++++++++++++

cd ../7up1
vi LTR7up1.dico
#HERVH-int    LTR7up1
grep "HERVH-int\|LTR7up1" ../../hg38.fa_notalt.out > LTR7up1.out

# put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7up1.out 



~/one_code_to_find_them_all.pl -rm LTR7up1.out -ltr LTR7up1.dico > LTR7up1.OC



cat *ltr.csv > LTR7up1.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7up1.csv -ltr LTR7up1=455 -int HERVH-int -ilen 7713 -rn LTR7up1
grep '2ltr_int' LTR7up1.csv.newnames.cordinates.txt > LTR7up1_full.bed
grep '1ltr_int' LTR7up1.csv.newnames.cordinates.txt > LTR7up1_truncated.bed
grep 'Sltr' LTR7up1.csv.newnames.cordinates.txt > LTR7up1_solo.bed
wc -l *7up1_*.bed
  175 LTR7up1_full.bed
  117 LTR7up1_solo.bed
   86 LTR7up1_truncated.bed
  378 total
#560 via grep -c "LTR7up1" ../hg38.fa.out

#++++++++++++++++++++++++++++++++++++++++++++++

cd ../7up2
vi LTR7up2.dico
#HERVH-int    LTR7up2
grep "HERVH-int\|LTR7up2" ../../hg38.fa_notalt.out > LTR7up2.out

# put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7up2.out 



~/one_code_to_find_them_all.pl -rm LTR7up2.out -ltr LTR7up2.dico > LTR7up2.OC



cat *ltr.csv > LTR7up2.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7up2.csv -ltr LTR7up2=455 -int HERVH-int -ilen 7713 -rn LTR7up2
grep '2ltr_int' LTR7up2.csv.newnames.cordinates.txt > LTR7up2_full.bed
grep '1ltr_int' LTR7up2.csv.newnames.cordinates.txt > LTR7up2_truncated.bed
grep 'Sltr' LTR7up2.csv.newnames.cordinates.txt > LTR7up2_solo.bed
wc -l *7up2_*.bed
  19 LTR7up2_full.bed
  52 LTR7up2_solo.bed
  25 LTR7up2_truncated.bed
  96 total
#119 via grep -c "LTR7up2" ../hg38.fa.out

#++++++++++++++++++++++++++++++++++++++++++++++

cd ../7y
vi LTR7y.dico
#HERVH-int    LTR7YY_cons
grep "HERVH-int\|LTR7YY_cons" ../../hg38.fa_notalt.out > LTR7y.out

# put in erv info
sed -i 's/Unspecified/LTR\/ERV1/g' LTR7y.out 



~/one_code_to_find_them_all.pl -rm LTR7y.out -ltr LTR7y.dico > LTR7y.OC



cat *ltr.csv > LTR7y.csv
perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7y.csv -ltr LTR7YY_cons=439 -int HERVH-int -ilen 7713 -rn LTR7YY_cons
grep '2ltr_int' LTR7y.csv.newnames.cordinates.txt > LTR7y_full.bed
grep '1ltr_int' LTR7y.csv.newnames.cordinates.txt > LTR7y_truncated.bed
grep 'Sltr' LTR7y.csv.newnames.cordinates.txt > LTR7y_solo.bed
wc -l *7y_*.bed
  28 LTR7y_full.bed
  295 LTR7y_solo.bed
   50 LTR7y_truncated.bed
  373 total
#513 via grep -c "LTR7YY_cons" ../hg38.fa.out


