# Written 03/16/2021 by Jason David Chobirko

# This R script - LTR7_blockFilt.v2.r - will be used alongside the bash script LTR7_muscleCons.v2.sh. It 
# will input the muscle .afa output and then output X files for each subfamily which correlate with 
# the blocks given to you by Thomas based on the LTR7 subfamily consensus alignment. 
# Block borders (w/o Gaps):
#       1   2a  2b  3   4   5   6  
# 7b    84  127 167 258 308 398 
# 7c    102 167 207 258 318 405 
# 7y    84  127 167 258 319 406 
# 7o    61  110 X   154 215 302 
# 7bc   83  132 X   176 237 324 
# 7d1   94  145 X   199 261 348 
# 7d2   93  144 X   198 260 347 
# 7u2   84  136 176 220 281 368 
# 7u1   84  136 184 239 301 388 
# 7up2  84  136 183 238 300 387 
# 7up1  84  136 183 238 300 387 

# Load library for sequences
library(Biostrings); library(seqinr)

names <- c("LTR7B", "LTR7C", "7bc", "7o", "7d1", "7d2", "7u2", "LTR7Y", "7u1", "7up2", "7up1")

for (samp in 1:length(names)) {
  assign(paste("l", names[samp], "seq", sep = "_"), readDNAStringSet(paste(names[samp], "comb.afa", sep = "_")))
  seq_name <- names(eval(parse(text = paste("l", names[samp], "seq", sep = "_"))))
  sequence <- paste(eval(parse(text = paste("l", names[samp], "seq", sep = "_"))))
  assign(paste("align", names[samp], "seq", sep = "_"), t(as.data.frame(list(strsplit(as.character(data.frame(seq_name, sequence)$sequence), split = "")))))
  # row.names(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))) <- seq_name
  
  # Determine which row is the one with the consensus sequence
  conRow <- which(seq_name == names[samp])
  
  # Now initiate the main method to create each of the files as needed.
  count <- 0; curPos <- 1
  
  # Initialize the right break points for each names variable
  if (names[samp] == "LTR7B") {brkPt <- c(84, 127, 167, 258, 308, 398, sum(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, ] != "-") - 1) }
  if (names[samp] == "LTR7C") {brkPt <- c(102, 167, 207, 258, 318, 405, sum(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, ] != "-") - 1) }
  if (names[samp] == "LTR7Y") {brkPt <- c(84, 127, 167, 258, 319, 406, sum(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, ] != "-") - 1) }
  if (names[samp] == "7o") {brkPt <- c(61, 110, 110, 154, 215, 302, sum(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, ] != "-") - 1) }
  if (names[samp] == "7bc") {brkPt <- c(83, 132, 132, 176, 237, 324, sum(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, ] != "-") - 1) }
  if (names[samp] == "7d1") {brkPt <- c(94, 145, 145, 199, 261, 348, sum(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, ] != "-") - 1) }
  if (names[samp] == "7d2") {brkPt <- c(93, 144, 144, 198, 260, 347, sum(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, ] != "-") - 1) }
  if (names[samp] == "7u2") {brkPt <- c(84, 136, 176, 220, 281, 368, sum(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, ] != "-") - 1) }
  if (names[samp] == "7u1") {brkPt <- c(84, 136, 184, 239, 301, 388, sum(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, ] != "-") - 1) }
  if (names[samp] == "7up2") {brkPt <- c(84, 136, 183, 238, 300, 387, sum(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, ] != "-") - 1) }
  if (names[samp] == "7up1") {brkPt <- c(84, 136, 183, 238, 300, 387, sum(eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, ] != "-") - 1) }
  
  # Initialize the break point variables needed below. 
  prev_bp <- curBlock <-  0; cur_bp <- brkPt[1]; 
  
  while (count < brkPt[7]) {  # ncol(eval(parse(text = paste("align", names[samp], "seq", sep = "_"))))) {
    
    # Iterate along the row with the consensus sequence and only append count if the base is 
    # not a gap but append curPos either way. 
    if (eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, curPos] != "-") {
      count <- count + 1 
    } 
    
    # Regardless, move along the sequence. 
    curPos <- curPos + 1
    
    # Check to see if the curPos has reached one of the breakpoints!
    if (count %in% brkPt & eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[conRow, curPos] != "-") { 
      curBlock <- curBlock + 1
      
      if (curBlock == 3 & names[samp] %in% c("7o", "7bc", "7d1", "7d2")) {
        curBlock <- 4
      }
       
      # Initialize the name and sequence vectors!
      name_vec <- c(); seq_vec <- c()
      
      # You reached a subfamily-specific break point, perform the for loop to generate the 
      # vectors necessary to use seqinr to output a .fa file!
      for (row in 1:nrow(eval(parse(text = paste("align", names[samp], "seq", sep = "_"))))) {
        # First, check to see if the current row is the consensus row. If not, continue
        if (row != conRow) {
          # Then, see if the current row's block has any sequence in it. If yes, append the necessary 
          # data to sequence and name vector
          tmp <- paste(gsub('\\-', "", eval(parse(text = paste("align", names[samp], "seq", sep = "_")))[row, c(prev_bp:curPos)]), collapse = "")
          
          if (width(tmp) > 9) {
            
            name_vec <- c(name_vec, seq_name[row]); seq_vec <- c(seq_vec, tmp)
          } else {
            tmp <- ""
          }
        }
      }
      
      # You have reached the end of the file and should have all the block sequence and names! 
      # Use seqinr to write the .fa file but only if one vector has > 10 elements!
      if (length(name_vec) > 10) {
        write.fasta(as.list(seq_vec), name_vec, paste(names[samp], "_block", curBlock, ".preH.fa", sep = ""))
      }
      
      # Now make sure to change the previous break point and the next break point
      prev_bp <- curPos; cur_bp <- brkPt[curBlock]
    }
  }
}
