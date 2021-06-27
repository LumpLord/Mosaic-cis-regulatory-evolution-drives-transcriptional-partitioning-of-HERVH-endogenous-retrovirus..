#OneCodeToFindThemAll runs of ORIGINAL RepeatMasker version 4.0.5 repeat Library 20140131 definitions
#"dico" files like dico-homo_sapiens-20120124_LTR7b.txt are a .tab with 1 column containin "HERVH-int" and one with one HERVH LTR eg. "LTR7B"

cd /workdir/tc643/RMout/hg38_oneCode/7B

cp ../dico-homo_sapiens-20120124_LTR7b.txt .
cp /workdir/tc643/RMout/hg38.out .
grep "HERVH-int\|LTR7B" hg38.out > LTR7B.out
~/one_code_to_find_them_all.pl -rm LTR7B.out -ltr dico-homo_sapiens-20120124_LTR7b.txt > hg38_LTR7b.output
cat *ltr.csv > LTR7b.csv

perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7b.csv -ltr LTR7B=464 -int HERVH-int -ilen 7713 -rn LTR7B  
grep '2ltr_int' LTR7b.csv.newnames.cordinates.txt > LTR7B_full.bed
grep '1ltr_int' LTR7b.csv.newnames.cordinates.txt > LTR7B_truncated.bed
grep 'Sltr' LTR7b.csv.newnames.cordinates.txt > LTR7B_solo.bed
wc -l *7B*.bed 
113 LTR7B_full.bed
  524 LTR7B_solo.bed
   95 LTR7B_truncated.bed
  732 total


cd /workdir/tc643/RMout/hg38_oneCode/7C
cp ../dico-homo_sapiens-20120124_LTR7c.txt .
grep "HERVH-int\|LTR7C" /workdir/tc643/RMout/hg38.out > LTR7C.out
~/one_code_to_find_them_all.pl -rm LTR7C.out -ltr dico-homo_sapiens-20120124_LTR7c.txt > hg38_LTR7c.output
cat *ltr.csv > LTR7c.csv

perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7c.csv -ltr LTR7C=471 -int HERVH-int -ilen 7713 -rn LTR7C  
grep '2ltr_int' LTR7c.csv.newnames.cordinates.txt > LTR7C_full.bed
grep '1ltr_int' LTR7c.csv.newnames.cordinates.txt > LTR7C_truncated.bed
grep 'Sltr' LTR7c.csv.newnames.cordinates.txt > LTR7C_solo.bed
wc -l *7C*.bed 
24 LTR7C_full.bed
  223 LTR7C_solo.bed
   50 LTR7C_truncated.bed


cd /workdir/tc643/RMout/hg38_oneCode/7Y
cp ../dico-homo_sapiens-20120124_LTR7y.txt .
grep "HERVH-int\|LTR7Y" /workdir/tc643/RMout/hg38.out > LTR7Y.out
~/one_code_to_find_them_all.pl -rm LTR7Y.out -ltr dico-homo_sapiens-20120124_LTR7y.txt > hg38_LTR7y.output
cat *ltr.csv > LTR7y.csv

perl ~/rename_mergedLTRelements_v6.0.pl -f LTR7y.csv -ltr LTR7Y=472 -int HERVH-int -ilen 7713 -rn LTR7Y  
grep '2ltr_int' LTR7y.csv.newnames.cordinates.txt > LTR7Y_full.bed
grep '1ltr_int' LTR7y.csv.newnames.cordinates.txt > LTR7Y_truncated.bed
grep 'Sltr' LTR7y.csv.newnames.cordinates.txt > LTR7Y_solo.bed
wc -l *7Y*.bed 
77 LTR7Y_full.bed
   77 LTR7Y_solo.bed
   24 LTR7Y_truncated.bed