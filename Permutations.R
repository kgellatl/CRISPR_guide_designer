library(Biostrings)

IUPAC_CODE_MAP #note the dictionary of definitions
seq <- c("A", "N", "C", "C", "T") #define vector of bases

seq_iupac <- as.vector(IUPAC_CODE_MAP[seq]) #gets all possible values for the iupac codes

seqs <- c() #vector will be generated with all values

for (i in 1:length(seq_iupac)) {
  if( i == 1) {
    seqs <- unlist(strsplit(seq_iupac[i], split = NULL)) # first entry will not need appending
  } else {
    seqs_new <- expand.grid(seqs, unlist(strsplit(seq_iupac[i], split = NULL))) # now to that intial vector we add all permutations
    seqs_new <- apply(seqs_new,1,paste0, collapse = "") #collapse the result into strings
    seqs <- seqs_new #update the seqs vector
  }
}

print(seqs)
print(unique(seqs))
