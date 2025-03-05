# ZAPbind
Code to analyse and recode CpG-containing sequences and optimal ZAP-binding sites in nucleic acid sequences.

## Overview
Collection of R scripts to analyse and recode DNA sequences. ZAP_Predicts.R calcualtes ZAP sensitivity scores of input sequences. RemovesCG.R synonymously recodes a DNA sequence by removing CpG dinucleotides while maintaining protein coding sequence. Recodes_to_EquiMono.R synonymously recodes a DNA sequence to contain a balance mononucleotides distribution (target 25% each).

## Prerequisites
This script requires R and the seqinr package.
To install seqinr, run:

```R
install.packages("seqinr")
```

## Usage
### Prepare input sequence
The script expects a FASTA file with Ensembl headers. Make sure *gene_symbol* information is included.
Example header format:
```ruby
>ENST00000421495.6 cds chromosome:GRCh38:1:1311606:1324637:-1 gene:ENSG00000127054.22 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:INTS11 description:integrator complex subunit 11 [Source:HGNC Symbol;Acc:HGNC:26052]
```
The gene_symbol: field (e.g., INTS11) is extracted.
## Run the script
Execute the script in R.

## Example Output
A snippet of ENST00000421495_6.csv:
Sequence_Name	Gene_Symbol	Length	CG_Number	CG_Frequency	CG_score	Distance_Raw	Distance_intra	Distance_score	Composition_score	ZS_Score
ENST00000421495.6	INTS11	1500	34	0.0227	0.1312	15	0.45	0.0075	0.62	0.0006
