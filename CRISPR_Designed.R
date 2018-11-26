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
  
  #Forward
  guides_f <- vector(mode = "list", length = nrow(nuclease_table))
  for (i in 1:nrow(nuclease_table)) {
    int_nuc <- nuclease_table[i,]
    pam <- as.character(int_nuc$PAM)
    pam <- DNAString(pam)
    for (j in 1:ncol(dat)) {
      int_seq <- dat[,j]
      int_seq <- paste0(int_seq, collapse = "")
      int_seq <- DNAString(int_seq)
      occurances <- matchPattern(pam, int_seq, fixed = FALSE)
      start_pos <- start(occurances@ranges)
      start_pos_cut <- start_pos + int_nuc$Cut # FIND CUT POSITION
      for (k in 1:length(mut_pos)) {
        int_mut_pos <- mut_pos[k]
        start_pos_mut <- start_pos_cut - int_mut_pos # DISTANCE OF CUT POSITION FROM MUTANT POSITIONS
        good <- which(abs(start_pos_mut)<=stringency) # WHICH ARE CLOSER THAN STRINGENCY
        if(length(good >= 1)){
          for (l in 1:length(good)) {
            goodint <- good[l]
            cut_dist <- start_pos_mut[good][l]
            if(cut_dist > 0){
              cut_dist <- cut_dist-1
            }
            int_new <- occurances[goodint]
            if(int_nuc$Cut < 0){
              if(start(int_new)-int_nuc$Guide.Length > 0){
                if(end(int_new) < nchar(int_seq)){
                  if( start(occurances)[goodint] <= int_mut_pos &&  int_mut_pos <= end(occurances)[goodint]){
                    dist <- 0
                  } else {
                    dist <-  start_pos[good][l]-int_mut_pos
                  }
                  guide <- subseq(int_seq, start = start(int_new)-int_nuc$Guide.Length, end = end(int_new)-nchar(pam))
                  guide_pam <- subseq(int_seq, start = start(int_new), end = end(int_new))
                  guide <- as.character(guide)
                  guide <- paste0(guide, "_", colnames(dat)[j], "_", guide_pam, "_", cut_dist, "_", dist)
                  guides_f[[i]] <- append(guides_f[[i]], guide)
                }
              }
            } else {
              if(start(int_new) > 0){
                if(end(int_new)+int_nuc$Guide.Length < nchar(int_seq)){
                  if( start(occurances)[goodint] <= int_mut_pos &&  int_mut_pos <= end(occurances)[goodint]){
                    dist <- 0
                  } else {
                    dist <-  int_mut_pos-start_pos[good][l]
                  }
                  guide <- subseq(int_seq, start = start(int_new)+nchar(pam), end = end(int_new)+int_nuc$Guide.Length)
                  guide_pam <- subseq(int_seq, start = start(int_new), end = end(int_new))
                  guide <- as.character(guide)
                  guide <- paste0(guide, "_", colnames(dat)[j], "_", guide_pam, "_", cut_dist, "_", dist)
                  guides_f[[i]] <- append(guides_f[[i]], guide)
                }
              }
            }
          }
        }
      }
    }
  }
  
  #Reverse
  guides_r <- vector(mode = "list", length = nrow(nuclease_table))
  for (i in 1:nrow(nuclease_table)) {
    int_nuc <- nuclease_table[i,]
    pam <- as.character(int_nuc$PAM)
    pam <- DNAString(pam)
    for (j in 1:ncol(dat)) {
      int_seq <- dat[,j]
      int_seq <- paste0(int_seq, collapse = "")
      int_seq <- DNAString(int_seq)
      int_seq <- reverseComplement(int_seq)
      occurances <- matchPattern(pam, int_seq, fixed = FALSE)
      start_pos <- start(occurances@ranges)
      start_pos_cut <- start_pos + int_nuc$Cut # FIND CUT POSITION
      for (k in 1:length(mut_pos)) {
        int_mut_pos <- mut_pos[k]
        int_mut_pos <- length(dat[,1]) - int_mut_pos +1
        start_pos_mut <- start_pos_cut - int_mut_pos # DISTANCE OF CUT POSITION FROM MUTANT POSITIONS
        good <- which(abs(start_pos_mut)<=stringency) # WHICH ARE CLOSER THAN STRINGENCY
        if(length(good >= 1)){
          for (l in 1:length(good)) {
            goodint <- good[l]
            cut_dist <- start_pos_mut[good][l]
            cut_dist <- cut_dist-1
            int_new <- occurances[goodint]
            if(int_nuc$Cut < 0){
              if(start(int_new)-int_nuc$Guide.Length > 0){
                if(end(int_new) < nchar(int_seq)){
                  if( start(occurances)[goodint] <= int_mut_pos &&  int_mut_pos <= end(occurances)[goodint]){
                    dist <- 0
                  } else {
                    dist <-  start_pos[good][l]-int_mut_pos
                  }
                  guide <- subseq(int_seq, start = start(int_new)-int_nuc$Guide.Length, end = end(int_new)-nchar(pam))
                  guide_pam <- subseq(int_seq, start = start(int_new), end = end(int_new))
                  guide <- as.character(guide)
                  guide <- paste0(guide, "_", colnames(dat)[j], "_", guide_pam, "_", cut_dist, "_", dist)
                  guides_r[[i]] <- append(guides_r[[i]], guide)
                }
              }
            } else {
              if(start(int_new) > 0){
                if(end(int_new)+int_nuc$Guide.Length < nchar(int_seq)){
                  if( start(occurances)[goodint] <= int_mut_pos &&  int_mut_pos <= end(occurances)[goodint]){
                    dist <- 0
                  } else {
                    dist <-  int_mut_pos-start_pos[good][l]
                  }
                  guide <- subseq(int_seq, start = start(int_new)+nchar(pam), end = end(int_new)+int_nuc$Guide.Length)
                  guide_pam <- subseq(int_seq, start = start(int_new), end = end(int_new))
                  guide <- as.character(guide)
                  guide <- paste0(guide, "_", colnames(dat)[j], "_", guide_pam, "_", cut_dist, "_", dist)
                  guides_r[[i]] <- append(guides_r[[i]], guide)
                }
              }
            }
          }
        }
      }
    }
  }
  
  full_guides <- mapply(c, guides_f, guides_r, SIMPLIFY=FALSE)
  names(full_guides) <- nuclease_table$Name
  guide_table <- matrix(nrow = length(unlist(full_guides)))
  guide_table[,1] <- unlist(full_guides)
  guide_table <- as.data.frame(guide_table)
  guide_table[,2] <- names(unlist(full_guides))
  guide_table[,3] <- NA
  guide_table[match(unlist(guides_f), as.character(guide_table[,1])),3] <- "+"
  guide_table[match(unlist(guides_r), as.character(guide_table[,1])),3] <- "-"
  
  guide_table[,4] <- matrix(unlist(strsplit(as.character(guide_table[,1]), split = "_")), ncol = 5, byrow = T)[,2]
  guide_table[,5] <- matrix(unlist(strsplit(as.character(guide_table[,1]), split = "_")), ncol = 5, byrow = T)[,3]
  guide_table[,6] <- matrix(unlist(strsplit(as.character(guide_table[,1]), split = "_")), ncol = 5, byrow = T)[,4]
  guide_table[,7] <- matrix(unlist(strsplit(as.character(guide_table[,1]), split = "_")), ncol = 5, byrow = T)[,5]
  
  for (i in 1:nrow(guide_table)) {
    for (j in 1:length(unique(nuclease_table)$Name)) {
      int <- guide_table[i,"V2"]
      int2 <- unique(nuclease_table)$Name[j]
      int2 <- as.character(int2)
      ind <- grep(int2, int)
      if(length(ind) ==1){
        guide_table[i,"V2"] <- int2
      }
    }
  }
  
  if(ncol(dat) == 1){
    trimmed_guides <- guide_table
    trimmed_guides[,1] <- matrix(unlist(strsplit(as.character(guide_table[,1]), split = "_")), ncol = 5, byrow = T)[,1]
    trimmed_guides <- trimmed_guides[,-4]
    unique_guides <- unique(trimmed_guides[,1])
    trimmed_guides2 <- trimmed_guides[match(unique_guides, trimmed_guides[,1]),]
    trimmed_guides3 <- trimmed_guides2[,c(1,4,2,3)]
    colnames(trimmed_guides3) <- c("Spacer_seq", "PAM", "Nuclease", "Strand")
    trimmed_guides <- trimmed_guides3
    rownames(trimmed_guides) <- seq(1:nrow(trimmed_guides))
  } else {
    trimmed_guides <- guide_table
    trimmed_guides[,1] <- matrix(unlist(strsplit(as.character(guide_table[,1]), split = "_")), ncol = 5, byrow = T)[,1]
    trimmed_guides[,8] <- paste0(trimmed_guides[,1], trimmed_guides[,5])
    unique_guides <- unique(trimmed_guides[,8])
    trimmed_guides2 <- trimmed_guides[match(unique_guides, trimmed_guides[,8]),]
    for (i in 1:nrow(trimmed_guides2)){
      int <- trimmed_guides2[i,"V8"]
      ind <- grep(paste0("^", int, "$"), trimmed_guides$V8)
      name <- paste0(trimmed_guides$V4[ind], collapse = "_")
      trimmed_guides2$V4[i] <- name
    }
    trimmed_guides <- trimmed_guides2
    trimmed_guides <- trimmed_guides[,c(1,5,2,7,6,3,4)]
    colnames(trimmed_guides) <- c("Spacer_seq", "PAM", "Nuclease", "SNP_PAM", "Cut_SNP", "Strand", "Genotype")
    rownames(trimmed_guides) <- seq(1:nrow(trimmed_guides))
  }
  if(remove_wt == T){
    remov_ind <- grep(paste0("^", "WT", "$"), trimmed_guides$Genotype)
    trimmed_guides <- trimmed_guides[-remov_ind,]
  }
  return(trimmed_guides)
}

guides <- find_guides(input = IVS_II_654, nuclease_table = CRISPR_Nuclease_Table, stringency = 6)
write.csv(guides, file = "guides.csv")
