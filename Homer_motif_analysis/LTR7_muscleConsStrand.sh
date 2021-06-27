#!/bin/bash

# Writen 03/22/2021 by Jason David Chobirko

# This script - LTR7_muscleConsStrand.sh - will take a bunch of named files in the naming scheme "subfamily"_cons.fa and "subfamily".fa to 
# perform a dual cat and muscle alignment on the same subfamily file. Then it will split the files in R and then output the files to 
# perform an all vs all motif enrichment analysis on individual blocks of sequence between families. Then it will take those output files 
# and make pretty plots using ggplot2 in R. 

# Set the path for muscle
export PATH=/programs/muscle:$PATH

# Set the path for homer commands
export PATH=/workdir/$USER/homer/bin:$PATH

# Set up a vector with the names to automate things!
names=("LTR7B" "LTR7C" "7bc" "7o" "7d1" "7d2" "7u2" "LTR7Y" "7u1" "7up2" "7up1")

# Combine each list of sequences with the respective consensus and then align them and make sure they're all one line per sequence
# plus header
for NAME in "${names[@]}"
do
	cat /workdir/jdc397/0_orthTEs/1_ChIP/6_LTR7_tommy/0_allLTR7/"${NAME}"_cons.fa /workdir/jdc397/0_orthTEs/1_ChIP/6_LTR7_tommy/0_allLTR7/"${NAME}".fa | muscle -in - -out - -quiet | perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' > "${NAME}"_comb.afa &
done

# Signify that the script is currently waiting
echo "Currently waiting..."

# Use the 'wait' command to make sure all muscle runs are finished before moving onward to the R scripts!
wait 

# Signify that the script is no longer waiting!
echo "Starting Rscript..."

# You can use the command Rscript to automatically run an R script in Bash! This script will generate all the sequence 
# blocks for each family based on the right-most recombination block positions
Rscript LTR7_blockFilt.v2.r

# Signify that the R script is done!
echo "Rscript finished!"

# Now that the Rscript is complete, perform the homer all vs all analysis by taking the code from LTR7_nestedHomer2.sh!

# Signify that Homer2 work is beginning!
echo "Starting all vs all pairwise homer2 find..."

for ((l=1; l<3; l++)); do

	for ((k=1; k<8; k++)); do

		echo "Working on block "${k}"..."

		# 'flist' contains the name of each file that will be analyzed via pairwise homer2 known
		flist=(*"${k}".preH.fa)

		# Set up the first homer2 known output to use as a base for both if loops. 
		# Remember things are 0-indexed!

		# For + strand
		if [[ ($l == "1") ]]; then

			homer2 known -i "${flist[0]}" -b "${flist[0]}" -p 4 -strand "+" -m /workdir/jdc397/0_orthTEs/1_ChIP/3_hESC_H9_barakat2018/2_peaksAnalyzed/homer.motif -p 5 | sed '1d' | sort -k1,1 | cut -f1-2,6-7 | awk -v a="${flist[0]}" '{OFS="\t"; print $1, $2, a, $3, $4}' > tmp1p.txt
			homer2 known -i "${flist[0]}" -b "${flist[1]}" -p 4 -strand "+" -m /workdir/jdc397/0_orthTEs/1_ChIP/3_hESC_H9_barakat2018/2_peaksAnalyzed/homer.motif -p 5 | sed '1d' | awk -v a="${flist[0]}" -v b="${flist[1]}" '{OFS="\t"; print $0, a "|vs|" b}' > tmpA1p.txt
			
		# For - strand
		else

			homer2 known -i "${flist[0]}" -b "${flist[0]}" -p 4 -strand "-" -m /workdir/jdc397/0_orthTEs/1_ChIP/3_hESC_H9_barakat2018/2_peaksAnalyzed/homer.motif -p 5 | sed '1d' | sort -k1,1 | cut -f1-2,6-7 | awk -v a="${flist[0]}" '{OFS="\t"; print $1, $2, a, $3, $4}' > tmp1m.txt
			homer2 known -i "${flist[0]}" -b "${flist[1]}" -p 4 -strand "-" -m /workdir/jdc397/0_orthTEs/1_ChIP/3_hESC_H9_barakat2018/2_peaksAnalyzed/homer.motif -p 5 | sed '1d' | awk -v a="${flist[0]}" -v b="${flist[1]}" '{OFS="\t"; print $0, a "|vs|" b}' > tmpA1m.txt
			
		fi

		max="${#flist[@]}"
		for ((i=0; i<max; i++)); do
			for ((j=0; j<max; j++)); do
				file1="${flist[i]}"
				file2="${flist[j]}"
				
				if [[ ( "$file1" == "${flist[0]}" )  && ( "$file2" != "${flist[0]}" ) ]]; then
					
					if [[ ($l == "1") ]]; then

						homer2 known -i "$file1" -b "$file2" -strand "+" -m /workdir/jdc397/0_orthTEs/1_ChIP/3_hESC_H9_barakat2018/2_peaksAnalyzed/homer.motif -p 5 | sed '1d' | sort -k1,1 | cut -f8-9 | awk -v a="$file2" '{OFS="\t"; print a, $0}' > tmp2p.txt
						paste tmp1p.txt tmp2p.txt > tmp3p.txt
						cat tmp3p.txt > tmp1p.txt
						rm tmp3p.txt
					# For - strand
					else

						homer2 known -i "$file1" -b "$file2" -strand "-" -m /workdir/jdc397/0_orthTEs/1_ChIP/3_hESC_H9_barakat2018/2_peaksAnalyzed/homer.motif -p 5 | sed '1d' | sort -k1,1 | cut -f8-9 | awk -v a="$file2" '{OFS="\t"; print a, $0}' > tmp2m.txt
						paste tmp1m.txt tmp2m.txt > tmp3m.txt
						cat tmp3m.txt > tmp1m.txt
						rm tmp3m.txt
					fi
				fi

				echo "Files: $file1 $file2"

				if [[ ( "$file1" == "${flist[0]}" )  && ( "$file2" == "${flist[1]}" ) ]]; then
					echo "Files are: $file1 $file2 at Nothing"

				else
				
					if [[ ($l == "1") ]]; then
 
						homer2 known -i "$file1" -b "$file2" -strand "+" -m /workdir/jdc397/0_orthTEs/1_ChIP/3_hESC_H9_barakat2018/2_peaksAnalyzed/homer.motif -p 5 | sed '1d' | awk -v a="$file1" -v b="$file2" '{OFS="\t"; print $0, a "|vs|" b}' > tmpA2p.txt
						cat tmpA1p.txt tmpA2p.txt >> tmpA3p.txt
						cat tmpA3p.txt > tmpA1p.txt
						rm tmpA3p.txt
					# For - strand
					else

						homer2 known -i "$file1" -b "$file2" -strand "-" -m /workdir/jdc397/0_orthTEs/1_ChIP/3_hESC_H9_barakat2018/2_peaksAnalyzed/homer.motif -p 5 | sed '1d' | awk -v a="$file1" -v b="$file2" '{OFS="\t"; print $0, a "|vs|" b}' > tmpA2m.txt
						cat tmpA1m.txt tmpA2m.txt >> tmpA3m.txt
						cat tmpA3m.txt > tmpA1m.txt
						rm tmpA3m.txt
					fi
				fi
			done
		done

		
		if [[ ($l == "1") ]]; then

			# Now remove and finalize the files
			cat tmp1p.txt > block"${k}"_allFam_plusValues.txt
			rm tmp1p.txt tmp2p.txt

			sort -k3,3g tmpA1p.txt > block"${k}"_allEnrich_plusValues.txt
			rm tmpA1p.txt tmpA2p.txt
		# For - strand
		else

			# Now remove and finalize the files
			cat tmp1m.txt > block"${k}"_allFam_minusValues.txt
			rm tmp1m.txt tmp2m.txt

			sort -k3,3g tmpA1m.txt > block"${k}"_allEnrich_minusValues.txt
			rm tmpA1m.txt tmpA2m.txt
		fi
	done
done


# Say that you're making beautiful plots
echo "Making Beautiful Plots... "

# You can use the command Rscript to automatically run an R script in Bash! 
Rscript LTR7_blockFiltPlotsStrand.r

echo "Done!"

