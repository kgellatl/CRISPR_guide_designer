library(Biostrings)

find_guides <- function(input, nuclease_table, stringency = 5, remove_wt = T){
  input2 <- as.character(unlist(input$V2)) # I hate factors... should just take your input and turn it into a vector
  dat <- matrix(ncol = length(input2), nrow = nchar(input2)[1]) # Now we make a dataframe
  colnames(dat) <- input$V1
  # Below will turn your input sequences into a df, where each row is a base position, and each column is WT / Mut1 / Mut2...
  for (i in 1:length(input2)) {
    int <- input2[i]
    dat[,i] <- unlist(strsplit(int, split = ""))
  }
  
  # Now we query WHERE in our sequence a SNP occurs
  if(ncol(dat) > 1){
    positions <- apply(dat,1,unique)
    mut_pos <- c()
    for (i in 1:length(positions)) {
      int <- positions[[i]]
      if(length(int) > 1){
        mut_pos <- c(mut_pos, i)
      }
    }
  } else {
    mut_pos <- seq(1:nrow(dat))
  }
  
  for (i in 1:ncol(dat)) {
    int_seq <- dat[,i]
    int_seq <- paste0(int_seq, collapse = "")
    int_seq <- DNAString(int_seq)
    int_seq <- reverseComplement(int_seq)
    int_seq <- unlist(strsplit(as.character(int_seq), split = ""))
    dat <- cbind(dat, int_seq)
    colnames(dat)[ncol(dat)] <- paste0("RC_", colnames(dat)[i])
  }
  
  guides <- c()
  guides_f <- vector(mode = "list", length = nrow(nuclease_table))
  guides_r <- vector(mode = "list", length = nrow(nuclease_table))
  
  #####################################
  # ITERATE THROUGH EACH ROW OF THE NUCLEASE TABLE
  #####################################
  for (i in 1:nrow(nuclease_table)) {  
    int_nuc <- nuclease_table[i,]
    pam <- as.character(int_nuc$PAM)
    pam <- DNAString(pam)
    nuclease <- int_nuc$Name
    cut_from_pam <- int_nuc$Cut
    guide_size <- int_nuc$Guide.Length
    #####################################
    # ITERATE THROUGH EACH SEQUENCE
    #####################################
    for (j in 1:ncol(dat)) {  
      int_seq <- dat[,j]
      int_seq <- paste0(int_seq, collapse = "")
      int_seq <- DNAString(int_seq)
      occurances <- matchPattern(pam, int_seq, fixed = FALSE)
      #####################################
      # ITERATE THROUGH EACH PAM IF FOUND
      #####################################
      if(length(start(occurances@ranges)) > 0){
        for (k in 1:length(start(occurances@ranges))) {  
          int_pam <- occurances@ranges[k]
          start_pam <- start(int_pam)  # PAM START
          end_pam <- end(int_pam)  # PAM END
          #####################################
          # ITERATE THROUGH EACH MUTANT POSITION 
          #####################################
          for(l in 1:length(mut_pos)){ 
            if(j <= dim(input)[1]){  # MUTANT POSITION FOR + STRAND
              int_mut <- mut_pos[l] 
            } else {  # MUTANT POSITION FOR - STRAND
              int_mut <- (nrow(dat)-mut_pos[l])+1
            }
            
            #####################################
            # THIS IS FOR THE MINUS CUTTERS CAS9
            #####################################
            if(cut_from_pam < 0){
              cut_pos <- start_pam+cut_from_pam
              if(abs(cut_pos-int_mut) <= stringency){ #CHECK STRINGENCY
                if((start_pam-guide_size) > 0){ #CHECK GUIDE BOUNDARY
                  cut_dist <- cut_pos-int_mut
                  if(abs(cut_dist) == 1){
                    cut_dist <- 0
                  }
                  if(end_pam < int_mut){ #PAM TO LEFT OF MUT
                    pam_dist <- int_mut-end_pam
                  }
                  if(start_pam > int_mut){ #PAM TO RIGHT OF MUT
                    pam_dist <- (start_pam-int_mut)-1
                  }
                  if(mut_pos %in% seq(start_pam, end_pam)){ # MUT IN PAM
                    pam_dist <- 0
                  }
                  guide_start <- start_pam-guide_size
                  guide_end <- start_pam-1
                  seq <- subseq(int_seq, guide_start, guide_end)
                  seq <- as.character(seq)
                  pam_seq <- subseq(int_seq, start_pam, end_pam)
                  pam_seq <- as.character(pam_seq)
                  nuclease <- as.character(nuclease)
                  if(j <= dim(input)[1]){
                    strand <- "+"
                  } else {
                    strand <- "-"
                  }
                  genotype <- colnames(dat)[j]
                  guides <- c(guides, c(seq, pam_seq, nuclease, pam_dist, cut_dist, strand, genotype, paste(i,j,k,l)))
                }
              }
              #####################################
              # THIS IS FOR THE PLUS CUTTERS CPF1
              #####################################
            } else { 
              
              
              
              
            }
          }
        }
      }
    }
  }
  
  #####################################
  # BEGIN GUIDE PARSING
  #####################################
  guides <- as.data.frame(matrix(guides, ncol = 8, byrow = T))
  colnames(guides) <- c("Spacer", "PAM", "Nuclease", "SNP_PAM", "Cut_SNP", "Strand", "Genotype")
  minus <- which(guides$Strand == "-")
  guides$Genotype[minus] <- matrix(unlist(strsplit(as.character(guides$Genotype[minus]), split = "_")), ncol = 2, byrow = T)[,2]
  View(guides)
}
