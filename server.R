library(Biostrings)
library(Rfastp)
library(ggplot2)
library(data.table)
library(stringr)
library(tibble)
library(dplyr)
library(bslib)

print("WHOOOOOO")

merging <- function(file1, file2, out) {
  o <- rfastp(read1 = file1, read2 = file2,
              merge = TRUE, outputFastq = "", mergeOut = out, thread = 1)
}

convert_StringSet_to_df = function(stringset, combine_duplicates = T) {
  sequence_df = as_tibble(as.character(stringset))
  names(sequence_df) = "sequence"
  if (combine_duplicates) {
    sequence_df = sequence_df %>%
      group_by(sequence) %>%
      summarise(read_count = n()) %>%
      ungroup() %>%
      arrange(desc(read_count), sequence)
  } else {
    sequence_df = sequence_df %>%
      mutate(read_count = 1)
  }
  return(as.data.frame(sequence_df))
}

translate_reads = function(reads) {
  translations = translate(reads, if.fuzzy.codon = "solve")
  return(as.character(translations))
}

temp_output <- "./temp_fastq_output.fastq.gz"

merging("./r1_fastq.fastq.gz", "./r2_fastq.fastq.gz", temp_output)
    
fasta <- readLines(temp_output)
      
file.remove(temp_output)

r1 <- readLines("./r1_fastq.fastq.gz")
r2 <- readLines("./r2_fastq.fastq.gz")

head <- r1[seq(1, length(r1), by = 4)]
r1 <- r1[seq(2, length(r1), by = 4)]
r2 <- as.vector(reverseComplement(DNAStringSet(r2[seq(2, length(r2), by = 4)])))

hi <- data.table(id = fasta[seq(1, length(fasta), by = 4)], overlap = fasta[seq(2, length(fasta), by = 4)])
hi$r1 <- r1[which(head %in% hi$id)]
hi$r2 <- r2[which(head %in% hi$id)]

hi[, combined := paste0(hi$r1, substr(hi$r2, nchar(hi$overlap) + 1, nchar(hi$r2)))][]

sequences <- hi$combined

primers <- fread("./primes.csv", header=TRUE)

out <- NULL
for (row in 1:nrow(primers)) {
  F_seq <- as.character(primers[row,"Fwd Seq"])
  R_seq <- as.character(primers[row,"Rev Seq"])
  
  n1 <- unlist(str_extract_all(sequences, paste0("^", F_seq, "(.*)", as.vector(reverseComplement(DNAStringSet(R_seq)))[1], "$")))
  n2 <- DNAStringSet(as.vector(reverseComplement(DNAStringSet(unlist(str_extract_all(sequences, paste0("^", R_seq, "(.*)", as.vector(reverseComplement(DNAStringSet(F_seq)))[1], "$")))))))
  
  n1 <- substring(n1, nchar(F_seq) + 1, nchar(n1) - nchar(R_seq))
  n2 <- substring(n2, nchar(R_seq) + 1, nchar(n2) - nchar(F_seq))
  
  n1 <- n1[n1 != ""] %>% convert_StringSet_to_df()
  n2 <- n2[n2 != ""] %>% convert_StringSet_to_df()
  
  colnames(n1)[colnames(n1) == "read_count"] <- primers[row,]$Primer
  colnames(n2)[colnames(n2) == "read_count"] <- primers[row,]$Primer
  
  comb <- rbind(n1, n2)
  
  comb <- comb %>% group_by(sequence) %>% summarise(across(primers[row,]$Primer, ~sum(.x, na.rm = TRUE)))
  
  if (is.null(out)) {
    out <- comb
  } else {
    out <- merge(out, comb, all = TRUE, by = "sequence")
    out[is.na(out)] <- 0
  }
}

write.csv(out, file="dna_sequences.csv")

dna0 <- out
dna1 <- dna0
dna1$sequence <- substring(dna1$sequence, 2)
dna2 <- dna0
dna2$sequence <- substring(dna2$sequence, 3)

dna0$sequence <- translate_reads(DNAStringSet(dna0$sequence))
dna0$id <- 1:nrow(dna0)
dna0 <- dna0 %>% filter((str_ends(sequence, "TVSS") | str_ends(sequence, "TVSA")) & (str_starts(sequence, "Q") | str_starts(sequence, "E")))

dna1$sequence <- translate_reads(DNAStringSet(dna1$sequence))
dna1$id <- 1:nrow(dna1)
dna1 <- dna1 %>% filter((str_ends(sequence, "TVSS") | str_ends(sequence, "TVSA")) & (str_starts(sequence, "Q") | str_starts(sequence, "E")))

dna2$sequence <- translate_reads(DNAStringSet(dna2$sequence))
dna2$id <- 1:nrow(dna2)
dna2 <- dna2 %>% filter((str_ends(sequence, "TVSS") | str_ends(sequence, "TVSA")) & (str_starts(sequence, "Q") | str_starts(sequence, "E")))

dna <- bind_rows(dna0, dna1, dna2)

dna <- dna %>% distinct(id, .keep_all = TRUE)

dna$id <- NULL

dna <- dna %>% group_by(sequence) %>% summarise(across(where(~ !is.character(.) & !is.factor(.)), ~sum(.x, na.rm = TRUE)))

dna <- dna[dna$sequence != "", ]

write.csv(dna, file="./amino_sequences.csv")
