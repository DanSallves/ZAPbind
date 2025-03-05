## This script takes as input a fasta file containing a  sequence in the DNA format and calculates a ZAP-sensitivity score for that sequence

library(seqinr)

#Enter path to fasta file containing the DNA sequences
dna2 <- read.fasta("", as.string = FALSE)

seqs <- length(dna2)
ZAP_scores <- data.frame(matrix(ncol = 11, nrow = seqs)) 
colnames(ZAP_scores) <- c('Sequence_Name', 'Gene_Symbol', 'Length', 'CG_Number', 'CG_Frequency', 
                          'CG_score', 'Distance_Raw', 'Distance_intra', 
                          'Distance_score', 'Composition_score', 'ZS_Score')

#Main loop
for (k in 1:seqs) {
  
  # Sequence info
  header <- attributes(dna2[[k]])$Annot  
  header <- gsub("\n", " ", header)  
  seq_name <- strsplit(header, " ")[[1]][1]  
  ZAP_scores[k, "Sequence_Name"] <- names(dna2)[k]

  gene_symbol_match <- regmatches(header, regexpr("gene_symbol:[^ ]+", header))
  if (length(gene_symbol_match) > 0) {
    gene_symbol <- sub("gene_symbol:", "", gene_symbol_match)  
  } else {
    gene_symbol <- NA 
  }
  ZAP_scores[k, "Gene_Symbol"] <- gene_symbol 
  
  #Extract sequence - Make sure it is in the DNA format
  seq_char <- unlist(dna2[[k]])  
  mono <- table(seq_char)  
  di <- table(paste0(seq_char[-length(seq_char)], seq_char[-1]))  
  len <- length(seq_char)  
  
  #Calculates CG score
  ZAP_scores[k, 'Length'] <- len
  ZAP_scores[k, 'CG_Number'] <- ifelse("cg" %in% names(di), di["cg"], 0)  
  ZAP_scores[k, 'CG_Frequency'] <- ifelse("cg" %in% names(di), di["cg"] / len, 0)  
  CG_score <- ZAP_scores[k, 'CG_Frequency'] / 0.173
  ZAP_scores[k, 'CG_score'] <- CG_score
  
  if (!("cg" %in% names(di)) || di["cg"] < 2) {
    #Rejects if there are fewer than 2 CpG sites in the sequence
    ZAP_scores[k, 'Distance_Raw'] <- 0
    ZAP_scores[k, 'Distance_intra'] <- 0
    ZAP_scores[k, 'Distance_score'] <- 0
    ZAP_scores[k, 'Composition_score'] <- 0
    ZAP_scores[k, 'ZS_Score'] <- 0
  } else {
    #Calculates CpG positions and distances for other cases
    positions <- which(paste0(seq_char[-length(seq_char)], seq_char[-1]) == "cg")
    distances <- diff(positions)
    
    #Calculates Distance scores
    ZAP_scores[k, 'Distance_Raw'] <- sum(distances > 11 & distances < 33)
    if (length(distances) > 0) {
      ZAP_scores[k, 'Distance_intra'] <- sum(distances > 11 & distances < 33) / length(distances)
    } else {
      ZAP_scores[k, 'Distance_intra'] <- 0
    }
    
    ZAP_scores[k, 'Distance_score'] <- ZAP_scores[k, 'Distance_Raw'] / 2000
    
    #Calculates Composition scores
    if (length(distances) > 0) {
      compositions <- sapply(1:(length(positions) - 1), function(i) {
        frag <- seq_char[(positions[i] + 2):(positions[i + 1] - 1)]
        if (length(frag) > 0) {  
          frag_table <- table(frag)  
          return((frag_table["a"] + frag_table["t"]) / sum(frag_table))  
        }
        return(0)  
      })
      
      composition_score <- mean(compositions[compositions > 0], na.rm = TRUE)
    } else {
      composition_score <- 0
    }
    
    ZAP_scores[k, 'Composition_score'] <- composition_score
    
    #Calculates ZS_Score
    ZAP_scores[k, 'ZS_Score'] <- CG_score * ZAP_scores[k, 'Distance_score'] * composition_score
  }
}

#Exports as csv file
write.csv(ZAP_scores, "./Human_final_v5.csv", row.names = FALSE)
