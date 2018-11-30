library(Biostrings)

find_guides <- function(input, nuclease_table, stringency = 7, remove_wt = T, diagnostic = F){
  #####################################
  # Parse Input
  #####################################
  input2 <- as.character(unlist(input$V2))
  dat <- matrix(ncol = length(input2), nrow = nchar(input2)[1])
  colnames(dat) <- input$V1
  for (i in 1:length(input2)) {
    int <- input2[i]
    dat[,i] <- unlist(strsplit(int, split = ""))
  }
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
  ##### GET REVERSE COMPLEMENT AS WELL
  for (i in 1:ncol(dat)) {
    int_seq <- dat[,i]
    int_seq <- paste0(int_seq, collapse = "")
    int_seq <- DNAString(int_seq)
    int_seq <- reverseComplement(int_seq)
    int_seq <- unlist(strsplit(as.character(int_seq), split = ""))
    dat <- cbind(dat, int_seq)
    colnames(dat)[ncol(dat)] <- paste0("RC_", colnames(dat)[i])
  }
  #####################################
  # ITERATE THROUGH EACH ROW OF THE NUCLEASE TABLE
  #####################################
  guides <- c()
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
          start_pam <- start(int_pam)
          end_pam <- end(int_pam)
          #####################################
          # ITERATE THROUGH EACH MUTANT POSITION
          #####################################
          for(l in 1:length(mut_pos)){
            ##### MUTANT POSITION FOR + STRAND
            if(j <= dim(input)[1]){
              int_mut <- mut_pos[l]
              ##### MUTANT POSITION FOR - STRAND
            } else {
              int_mut <- (nrow(dat)-mut_pos[l])+1
            }
            #####################################
            # THIS IS FOR THE MINUS CUTTERS CAS9
            #####################################
            if(cut_from_pam < 0){
              cut_pos <- start_pam+cut_from_pam
              ##### CHECK GUIDE BOUNDARY
              if(start_pam-guide_size > 0){
                ##### PAM TO LEFT OF MUT
                if(end_pam < int_mut){
                  pam_dist <- int_mut-end_pam
                  cut_dist <- cut_pos-int_mut
                  ##### DISTANCES FOR EDGE CUTS
                  if(abs(cut_dist) == 1){
                    cut_dist <- 0
                  }
                }
                ##### PAM TO RIGHT OF MUT
                if(start_pam > int_mut){
                  pam_dist <- int_mut-start_pam
                  cut_dist <- (int_mut-cut_pos)
                  if(cut_pos > int_mut){
                    cut_dist <- cut_dist + 1
                  }
                }
                ##### MUT IN PAM
                if(int_mut %in% seq(start_pam, end_pam)){
                  pam_dist <- 0
                  cut_dist <- cut_pos-int_mut
                }
                ##### SELECT ONLY THOSE FOR STRINGENCY
                if(abs(cut_dist) <= stringency){
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
              cut_pos1 <- end_pam+cut_from_pam
              cut_pos2 <- end_pam+cut_from_pam+4
              if(cut_pos1 > int_mut){
                cut_pos <- cut_pos1
              }
              if(cut_pos2 < int_mut){
                cut_pos <- cut_pos2
              }
              if(cut_pos1 <= int_mut && cut_pos2 >= int_mut){
                cut_pos <- int_mut
              }
              ##### CHECK GUIDE BOUNDARY
              if(end_pam+guide_size < nrow(dat)){
                ##### PAM TO LEFT OF MUT
                if(end_pam < int_mut){
                  pam_dist <- int_mut-end_pam
                  cut_dist <- cut_pos-int_mut+1
                  if(cut_pos > int_mut){
                    cut_dist <- cut_dist-1
                  }
                }
                ##### PAM TO RIGHT OF MUT # VERY UNLIKELY BUT PROBABLY WRONT
                if(start_pam > int_mut){
                  pam_dist <- start_pam-int_mut
                  cut_dist <- cut_pos-int_mut
                }
                ##### MUT IN PAM # VERY UNLIKELY BUT PROBABLY WRONT
                if(int_mut %in% seq(start_pam, end_pam)){
                  pam_dist <- 0
                  cut_dist <- cut_pos-int_mut
                }
                ##### DISTANCES FOR EDGE CUTS
                if(abs(cut_dist) == 1){
                  cut_dist <- 0
                }
                ##### SELECT ONLY THOSE FOR STRINGENCY
                if(abs(cut_dist) <= stringency){
                  guide_start <- end_pam+1
                  guide_end <- end_pam+guide_size
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
            }
          }
        }
      }
    }
  }
  #####################################
  # FORMAT OUTPUT
  #####################################
  guides_format <- as.data.frame(matrix(guides, ncol = 8, byrow = T))
  colnames(guides_format) <- c("Spacer", "PAM", "Nuclease", "PAM_to_SNP", "SNP_to_Cut", "Strand", "Genotype")
  minus <- which(guides_format$Strand == "-")
  guides_format$Genotype[minus] <- matrix(unlist(strsplit(as.character(guides_format$Genotype[minus]), split = "_")), ncol = 2, byrow = T)[,2]
  ##### FIND UNIQUE GUIDES FOR SINGLE INPUT
  if(ncol(dat) == 2){
    guides_format <- guides_format[,c(1,2,3,6)]
    if(diagnostic == FALSE){
      guides_format <- guides_format[,c(1,2,3,6)]
    } else {
      guides_format <- guides_format[,c(1,2,3,6,8)]
    }
    full_guide <- apply(guides_format[,1:2],1,paste0, collapse = "")
    ogs <- unique(full_guide)
    keep <- match(ogs, full_guide)
    guides_format <- guides_format[keep,]
    ##### FIND UNIQUE GUIDES FOR MULTIPLE INPUT
  } else {
    if(diagnostic == FALSE){
      guides_format <- guides_format[,1:7]
    } else {
      guides_format <- guides_format[,1:8]
    }
    full_guide <- apply(guides_format[,1:2],1,paste0, collapse = "")
    ogs <- unique(full_guide)
    guides_format$Genotype <- as.character(guides_format$Genotype)
    remove <- c()
    for (i in 1:length(ogs)) {
      int <- ogs[i]
      ind <- grep(paste0("^",int,"$"), full_guide)
      if(length(ind) > 1){
        new_title <- paste0(guides_format$Genotype[ind], collapse = "_")
        guides_format$Genotype[ind] <- new_title
        remove <- c(remove, ind[2])
      }
    }
    guides_format <- guides_format[-remove,]
    ##### REMOVE ONLY WT TARGETTING
    if(remove_wt){
      remove_ind <- grep(paste0("^", "WT", "$"), guides_format$Genotype)
      guides_format <- guides_format[-remove_ind,]
    }
  }
  #####################################
  # FORMAT OUTPUT
  #####################################
  rownames(guides_format) <- seq(1:nrow(guides_format))
  print(table(guides_format$Nuclease))
  return(guides_format)
}
