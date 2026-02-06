library(Biostrings)
library(rentrez)
library(stringr)

# Read your CODH FASTA file
codh_seqs <- readAAStringSet("../cladeB-CODH/CladeB.fasta")  # or readDNAStringSet if nucleotide

# Extract organism names from headers (text within [ ])
extract_organism <- function(header) {
  # Extract text between [ and ]
  org <- str_extract(header, "(?<=\\[)[^\\]]+(?=\\])")
  return(org)
}

# Get all organism names
organisms <- sapply(names(codh_seqs), extract_organism)
organisms <- unique(organisms)  # Remove duplicates
organisms <- organisms[!is.na(organisms)]  # Remove NAs if any

print(paste("Found", length(organisms), "unique organisms"))
print(head(organisms))

# Function to search and download 16S
download_16S <- function(species_name) {
  cat("Searching for:", species_name, "\n")
  
  search_term <- paste0(species_name, "[ORGN] AND (16S ribosomal RNA[Title] OR 16S rRNA[Title])")
  
  tryCatch({
    search <- entrez_search(db = "nucleotide", 
                            term = search_term,
                            retmax = 1)
    
    if (search$count > 0) {
      fetch <- entrez_fetch(db = "nucleotide", 
                           id = search$ids[1], 
                           rettype = "fasta")
      cat("  ✓ Found:", search$ids[1], "\n")
      return(fetch)
    } else {
      cat("  ✗ No 16S found\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("  ✗ Error:", e$message, "\n")
    return(NULL)
  })
}

# Download all 16S sequences (with rate limiting)
sequences_16S <- list()
failed <- c()

for (i in seq_along(organisms)) {
  org <- organisms[i]
  cat(paste0("[", i, "/", length(organisms), "] "))
  
  Sys.sleep(0.4)  # Be nice to NCBI servers (max ~3 requests/second)
  
  result <- download_16S(org)
  
  if (!is.null(result)) {
    sequences_16S[[org]] <- result
  } else {
    failed <- c(failed, org)
  }
}

# Save results - CORRECTED
if (length(sequences_16S) > 0) {
  # Write each sequence to file
  cat(paste(sequences_16S, collapse = ""), file = "16S_sequences.fasta")
  
  # OR better, write line by line:
  writeLines(unlist(sequences_16S), "16S_sequences.fasta")
  
  cat("\n✓ Successfully downloaded", length(sequences_16S), "16S sequences\n")
}

if (length(failed) > 0) {
  cat("\n✗ Failed:", length(failed), "organisms\n")
  write.table(data.frame(organism = names(failed)), 
              "failed_16S_downloads.txt", 
              row.names = FALSE, quote = FALSE, col.names = FALSE)
}

# Summary
cat("\nSummary:")
cat("\n  Total organisms:", length(organisms))
cat("\n  Successful:", length(sequences_16S))
cat("\n  Failed:", length(failed))
cat("\n  Success rate:", round(100*length(sequences_16S)/length(organisms), 1), "%\n")

# Save successful downloads
if (length(sequences_16S) > 0) {
  cat(sequences_16S, file = "16S_sequences.fasta", sep = "")
  cat("\nSuccessfully downloaded", length(sequences_16S), "16S sequences\n")
} else {
  cat("\nNo sequences downloaded!\n")
}

# Report failures
if (length(failed) > 0) {
  cat("\nFailed to find 16S for", length(failed), "organisms:\n")
  cat(paste(failed, collapse = "\n"), "\n")
  write.table(data.frame(organism = failed), 
              "failed_16S_downloads.txt", 
              row.names = FALSE, quote = FALSE)
}

# Create a mapping file for later use
if (length(sequences_16S) > 0) {
  mapping_16S <- data.frame(
    organism = names(sequences_16S),
    has_16S = TRUE
  )
  
  write.csv(mapping_16S, "16S_organism_mapping.csv", row.names = FALSE)
}


### also donwload rpoB for the same organisms
# Download RpoB sequences (commonly used)
download_control_protein <- function(species_name, protein = "rpoB") {
  search_term <- paste0(species_name, "[ORGN] AND ", protein, "[Gene]")
  
  search <- entrez_search(db = "protein", 
                          term = search_term,
                          retmax = 1)
  
  if (search$count > 0) {
    fetch <- entrez_fetch(db = "protein", 
                         id = search$ids[1], 
                         rettype = "fasta")
    return(fetch)
  }
  return(NULL)
}

# Download all RpoB sequences (with rate limiting)
sequences_rpoB <- list()
failed_rpoB <- c()
for (i in seq_along(organisms)) {
  org <- organisms[i]
  cat(paste0("[", i, "/", length(organisms), "] "))
  
  Sys.sleep(0.4)  # Be nice to NCBI servers
  
  result <- download_control_protein(org, protein = "rpoB")
  
  if (!is.null(result)) {
    sequences_rpoB[[org]] <- result
    cat("  ✓ Found RpoB\n")
  } else {
    failed_rpoB <- c(failed_rpoB, org)
    cat("  ✗ No RpoB found\n")
  }
}
# Save RpoB results
if (length(sequences_rpoB) > 0) {
  writeLines(unlist(sequences_rpoB), "rpoB_sequences.fasta")
  cat("\n✓ Successfully downloaded", length(sequences_rpoB), "RpoB sequences\n")
}

# Maps CODH accesion number to rpoB presence/absence

# Phrase to extract accession number from CODH headers
extract_accesion <- function(header){
  acc <- str_extract(header, "^[^ ]+")  # Extract text before first space
  return(acc)
}

codh_accessions <- sapply(names(codh_seqs), extract_accesion)
codh_organisms <- sapply(names(codh_seqs), extract_organism)
# make df of accesion and organism name
codh_df <- data.frame(
  codh_accession = codh_accessions,
  organism = codh_organisms
)

sequences_rpoB_fasta <- readAAStringSet("rpoB_sequences.fasta")
sequences_16S_fasta <- readDNAStringSet("16S_sequences.fasta")

rpoB_accessions <- sapply(names(sequences_rpoB_fasta), extract_accesion)
rpoB_organisms <- sapply(names(sequences_rpoB_fasta), extract_organism)

# make df of accesion and organism name
rpoB_df <- data.frame(
  rpoB_accession = rpoB_accessions,
  organism = rpoB_organisms
)

extract_organism2 <- function(header) {
  parts <- stringr::str_split_fixed(header, " ", 4)
  return(paste(parts[, 2], parts[, 3], sep = " ")) 
}

S_accessions <- sapply(names(sequences_16S_fasta), extract_accesion)
S_organisms <- sapply(names(sequences_16S_fasta), extract_organism2)
# make df of accesion and organism name
sequences_16S_df <- data.frame(
  X16S_accession = S_accessions,
  organism = S_organisms
)

# merge codh_df with rpoB_df and sequences_16S_df by organism name
merged_df <- merge(codh_df, rpoB_df, by = "organism", all.x = TRUE)
merged_df <- merge(merged_df, sequences_16S_df, by = "organism", all.x = TRUE)

# save as csv
write.csv(merged_df, "CODH_rpoB_16S_mapping.csv", row.names = FALSE)
