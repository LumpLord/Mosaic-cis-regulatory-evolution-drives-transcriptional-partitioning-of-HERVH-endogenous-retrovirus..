#!/bin/bash

# Writen 09/23/2021 by Jason David Chobirko

# This script - LTR7_muscleConsStrandFind.sh - will take a bunch of named files in the naming scheme "subfamily"_cons.fa and "subfamily".fa to 
# perform a dual cat and muscle alignment on the same subfamily file. Then it will split the files in R and then output the files to 
# perform homer2 find on the 2b cluster against the Sox2/3 motif to determine whether the motif is enriched in active elements 

# Set the path for muscle
export PATH=/programs/muscle:$PATH

# Set the path for homer commands
export PATH=/workdir/$USER/homer/bin:$PATH

# Set up a vector with the names to automate things!
# names=$(ls *_cons.fa | sed 's/_cons.fa//g' | sort | uniq)
names=("LTR7B" "LTR7C" "LTR7bc" "LTR7o" "LTR7d1" "LTR7d2" "LTR7u2" "LTR7Y" "LTR7u1" "LTR7up2" "LTR7up1")


# Signify that the alignments are currently happening
echo "Begining alignments..."

# Combine each list of sequences with the respective consensus and then align them and make sure they're all one line per sequence
# plus header
for NAME in "${names[@]}"
do
	cat "${NAME}"_cons.fa "${NAME}".fa | muscle -in - -out - -quiet | perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' > "${NAME}"_comb.afa &
done

# Use the 'wait' command to make sure all muscle runs are finished before moving onward to the R scripts!
wait 

# Signify that the script is no longer waiting!
echo "Starting Rscript..."

# You can use the command Rscript to automatically run an R script in Bash! This script will generate all the sequence 
# blocks for each family based on the right-most recombination block positions
Rscript LTR7_blockFilt.r

# Signify that the R script is done!
echo "Rscript finished!"

# Now that the Rscript is complete, perform the homer all vs all analysis by taking the code from LTR7_nestedHomer2.sh!

# Signify that Homer2 work is beginning!
echo "Starting block3 (2b) homer2 find..."

# Use the .motifs file that has Sox3 in it on block3 (2b). 
for f in $(ls *block3.preH.fa); do homer2 find -i $f -m sox3Motif.motifs -p 8 > "${f:0:-8}"_find_2b.txt; done

# Now combine them all together into a single file! 
cat *_find_2b.txt | cut -f1-3,5-6 > comb_block2b_find.txt

# Finally, generate a stats file from all of the output find files to see what there is to see for info
printf '%s\t%s\t%s\t%s\t%s\n' "subfamily" "numTotal" "numWithSox3" "numSOLO" "num5UTR" > 2b_sox3_findStats.txt
for f in "${names[@]}"; do printf '%s\t%s\t%s\t%s\t%s\n' $f $(expr `cat "${f}"_block3.preH.fa | wc -l` / 2) $(cat "${f}"_block3_find_2b.txt | wc -l) $(grep $f "${f}"_block3_find_2b.txt | grep "SOLO" | wc -l) $(grep $f "${f}"_block3_find_2b.txt | grep "5UTR" | wc -l) >> 2b_sox3_findStats.txt; done

echo "Done!"

