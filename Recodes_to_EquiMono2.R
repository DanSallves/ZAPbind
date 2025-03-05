### - This script recodes an input coding sequence to a sequence with similar base 
### - composition (as close to 25% as possible) while maintaining protein sequence


codon_table <- list(
  A = c("GCT", "GCC", "GCA", "GCG"),
  C = c("TGT", "TGC"),
  D = c("GAT", "GAC"),
  E = c("GAA", "GAG"),
  F = c("TTT", "TTC"),
  G = c("GGT", "GGC", "GGA", "GGG"),
  H = c("CAT", "CAC"),
  I = c("ATT", "ATC", "ATA"),
  K = c("AAA", "AAG"),
  L = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
  M = c("ATG"),
  N = c("AAT", "AAC"),
  P = c("CCT", "CCC", "CCA", "CCG"),
  Q = c("CAA", "CAG"),
  R = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
  S = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
  T = c("ACT", "ACC", "ACA", "ACG"),
  V = c("GTT", "GTC", "GTA", "GTG"),
  W = c("TGG"),
  Y = c("TAT", "TAC"),
  STOP = c("TAA", "TAG", "TGA")
)

nucleotide_composition <- function(sequence) {
  total_length <- nchar(sequence)
  counts <- table(strsplit(sequence, "")[[1]])
  composition <- sapply(c("A", "C", "G", "T"), function(nuc) counts[nuc] / total_length)
  composition[is.na(composition)] <- 0  
  return(composition)
}


recode_sequence <- function(input_sequence) {
  codons <- substring(input_sequence, seq(1, nchar(input_sequence), 3), seq(3, nchar(input_sequence), 3))

  recoded_sequence <- ""
  
  for (codon in codons) {
    amino_acid <- names(which(sapply(codon_table, function(codons) codon %in% codons)))
    
    available_codons <- codon_table[[amino_acid]]
    
    current_composition <- nucleotide_composition(recoded_sequence)
    
    best_codon <- available_codons[1]
    best_score <- Inf
    
    for (test_codon in available_codons) {
      test_sequence <- paste0(recoded_sequence, test_codon)
      test_composition <- nucleotide_composition(test_sequence)
      
      balance_score <- sum((test_composition - 0.25)^2)
      
      if (balance_score < best_score) {
        best_score <- balance_score
        best_codon <- test_codon
      }
    }
    
    recoded_sequence <- paste0(recoded_sequence, best_codon)
  }
  
  return(recoded_sequence)
}

# Usage: paste sequence directly here or load from workign directory
input_sequence <- "ATGGCTTGTGATGAATTTGGTCATATTAAATTAAATCCTCAATCTACTTGGGTTTATATGGCTTGTGATGAATTTGGTCATATTAAATTAAATCCTCAATCTACTTGGGTTTATATGGCTTGTGATGAATTTGGTCATATTAAATTAAATCCTCAATCTACTTGGGTTTAT"
recoded_sequence <- recode_sequence(input_sequence)
print(paste("Recoded Sequence:", recoded_sequence))

# After recoding, you can verify nucleotide composition with this command
composition <- nucleotide_composition(recoded_sequence)
print("Nucleotide Composition:")
print(composition)