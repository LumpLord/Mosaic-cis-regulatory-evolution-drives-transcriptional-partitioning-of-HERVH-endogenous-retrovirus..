# Written 05/20/2021 by Jason David Chobirko

# This R script - LTR7_blockFiltPlotsStrand.r - will be used alongside the bash script 
# LTR7_muscleConsStrand.sh and the other R script LTR7_blockFilt.v2.r and will produce beautiful 
# plots combining the strand information into a single plot for each family.

library(plyr); library(tidyverse); library(ggrepel)

# Generate all of the plus/minus strand info and then plot it for each family comparison! You'll want to combine 
# the strand info for each family as well and not do any filtering to enrich for motifs that are mostly unique 
# to each family (those motifs that have the most enriched comparisons will have a higher median value which is 
# nice to see, I guess?) 

combEnrichPlus <- combEnrichMinus <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), c("motifName", "consensus", "pValue", "log_pValue", "enrichment", "numA", "percentA", "numB", "percentB", "temp"))

## Set this to the total number of blocks. 
for (block in 1:7) {
	# Create a combined data frame of all 7 output block files. 
	combEnrichPlus <- rbind(combEnrichPlus, read.table(file = paste("block", block, "_allEnrich_plusValues.txt", sep = ""), col.names = c("motifName", "consensus", "pValue", "log_pValue", "enrichment", "numA", "percentA", "numB", "percentB", "temp")))
	combEnrichMinus <- rbind(combEnrichMinus, read.table(file = paste("block", block, "_allEnrich_minusValues.txt", sep = ""), col.names = c("motifName", "consensus", "pValue", "log_pValue", "enrichment", "numA", "percentA", "numB", "percentB", "temp")))
}

# Add a new column and then fill it with the block numbers.
combEnrichPlus$block <- combEnrichMinus$block <- ""                                      
combEnrichPlus[grepl("block1", combEnrichPlus$temp), 11] <- 1; combEnrichMinus[grepl("block1", combEnrichMinus$temp), 11] <- 1
combEnrichPlus[grepl("block2", combEnrichPlus$temp), 11] <- 2; combEnrichMinus[grepl("block2", combEnrichMinus$temp), 11] <- 2
combEnrichPlus[grepl("block3", combEnrichPlus$temp), 11] <- 3; combEnrichMinus[grepl("block3", combEnrichMinus$temp), 11] <- 3
combEnrichPlus[grepl("block4", combEnrichPlus$temp), 11] <- 4; combEnrichMinus[grepl("block4", combEnrichMinus$temp), 11] <- 4
combEnrichPlus[grepl("block5", combEnrichPlus$temp), 11] <- 5; combEnrichMinus[grepl("block5", combEnrichMinus$temp), 11] <- 5
combEnrichPlus[grepl("block6", combEnrichPlus$temp), 11] <- 6; combEnrichMinus[grepl("block6", combEnrichMinus$temp), 11] <- 6
combEnrichPlus[grepl("block7", combEnrichPlus$temp), 11] <- 7; combEnrichMinus[grepl("block7", combEnrichMinus$temp), 11] <- 7

# Now to tidy up the data a bit. 
combEnrichPlus$temp <- as.character(combEnrichPlus$temp); combEnrichMinus$temp <- as.character(combEnrichMinus$temp)
combEnrichPlus$temp <- gsub("\\_block[0-7].preH.fa", "", x = combEnrichPlus$temp); combEnrichMinus$temp <- gsub("\\_block[0-7].preH.fa", "", x = combEnrichMinus$temp)
combEnrichPlus <- combEnrichPlus %>% separate(temp, into = c("famA", "famB"), sep = "[|]vs[|]"); combEnrichMinus <- combEnrichMinus %>% separate(temp, into = c("famA", "famB"), sep = "[|]vs[|]")
combEnrichPlus <- combEnrichPlus %>% separate(motifName, into = c("short", "rest"), sep = "[(]"); combEnrichMinus <- combEnrichMinus %>% separate(motifName, into = c("short", "rest"), sep = "[(]")

# Below is the combined version of the names variable. Yeah!
if (length(unique(combEnrichPlus$famA)) < 10) {
	names <- c("comb_bco", "comb_d1_2", "7u1", "7u2", "comb_up1_2", "LTR7B", "LTR7C", "LTR7Y"); pdf(file = "052021_famMotEnrCombPlus.pdf")
	} else { 
	names <- c("7bc", "7d1", "7d2", "7o", "7u1", "7u2", "7up1", "7up2", "LTR7B", "LTR7C", "LTR7Y"); pdf(file = "052021_famMotEnrPlus.pdf") }

combEnrichPlus$block <- as.factor(combEnrichPlus$block); combEnrichMinus$block <- as.factor(combEnrichMinus$block)

# Perform the main plotting loop for each family!
for (num in 1:length(names)) {
	tempP <- combEnrichPlus %>% subset(famA == names[num]) %>% arrange(pValue)
	tempM <- combEnrichMinus %>% subset(famA == names[num]) %>% arrange(pValue)
 
	## This will iterate over each group and generate a summarized version of each motif 
	## enrichment found there for ease of plotting, hopefully?
	for (g in 1:length(levels(tempP$block))) {
	
		temp2P <- tempP %>% subset(block == levels(tempP$block)[g])
		temp2M <- tempM %>% subset(block == levels(tempM$block)[g])

		if (nrow(temp2P) == 0) {next}		

		## Generate the list of unique motifs present in the current group
		ntmpP <- unique(temp2P$short); ntmpM <- unique(temp2M$short)
			
		## Now iterate over the number of motifs found in the current positive group
		for (row in 1:length(ntmpP)) {
			  			  
			# Generate the data.frame to store all the median p-values for each family in! Yeah!
			if (g == 1 & row == 1) {
				tempFrameP <- data.frame(ntmpP[row], median(subset(temp2P, short == ntmpP[row])$log_pValue), levels(tempP$block)[g]); colnames(tempFrameP) <- c("ID", "log_pValue", "block") 
			} else {
				thing <- data.frame(ntmpP[row], median(subset(temp2P, short == ntmpP[row])$log_pValue), levels(tempP$block)[g]); colnames(thing) <- c("ID", "log_pValue", "block")
				tempFrameP <- rbind(tempFrameP, thing)
			}
		}

		## Now iterate over the number of motifs found in the current negative group
		for (row in 1:length(ntmpM)) {
			  			  
			# Generate the data.frame to store all the median p-values for each family in! Yeah!
			if (g == 1 & row == 1) {
				tempFrameM <- data.frame(ntmpM[row], median(subset(temp2M, short == ntmpM[row])$log_pValue), levels(tempM$block)[g]); colnames(tempFrameM) <- c("ID", "log_pValue", "block") 
			} else {
				thing <- data.frame(ntmpM[row], median(subset(temp2M, short == ntmpM[row])$log_pValue), levels(tempM$block)[g]); colnames(thing) <- c("ID", "log_pValue", "block")
				tempFrameM <- rbind(tempFrameM, thing)
			}
		}
	}

# Once both tables have been generated, rbind them together and make the plots!
tempFrameP$strand <- "+"; tempFrameM$strand <- "-"
comb <- rbind(tempFrameP, tempFrameM)

# Change the block labels to match the recombination naming scheme. 
comb$block <- revalue(comb$block, c("1" = "1", "2"="2a", "3"="2b", "4"="3", "5"="4", "6"="5", "7"="6"))
		
## Now save the final table and print it into the pdf!
print(ggplot(data = comb, aes(x = block, y = -log_pValue, color = strand)) + geom_point(alpha = 0.80) + guides(alpha = FALSE) + scale_color_brewer(palette = "Set1") + geom_text_repel(data = subset(comb, log_pValue <= quantile(comb$log_pValue, .005)), aes(x = block, y = -log_pValue, color = strand, label = ID), max.iter = 100, size = 3.25) + theme_bw() + labs(title = paste(names[num], "Family Motif Enrichment", sep = " "), x = "Consensus Sequence Block", y = "Median Motif Enrichment -log(pValue)") + theme(plot.title = element_text(size = 20, face = "bold.italic"), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14), axis.text.y = element_text(size = 14)))		

## Finally, output the sumTable!
comb <- comb %>% arrange(block, log_pValue)
write.table(comb, file = paste("comb_", names[num], ".txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
}
dev.off()
