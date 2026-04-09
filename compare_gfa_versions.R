# Comparison script: GFA C-program vs gfaR R package
library(gfaR)

# 1. Setup paths
fasta_file <- "test_data/gfa_test.fasta"
c_program <- "test_data/gfa.exe"
output_prefix <- "test_data/c_output"

# 2. Run the original C program
cat("--- Running original C program ---\n")
system(paste(c_program, "-seq", fasta_file, "-chrom Test -out", output_prefix))

# 3. Run the R package
cat("\n--- Running gfaR R package ---\n")
read_fasta <- function(file) {
    lines <- readLines(file)
    seq_lines <- lines[!grepl("^>", lines)]
    paste0(seq_lines, collapse = "")
}
dna_seq <- read_fasta(fasta_file)
r_results <- gfa_analyze(dna_seq)

# 4. Define comparison function
compare_results <- function(r_df, gff_file, motif_name) {
    cat(sprintf("\nComparing %s:\n", motif_name))
    
    if (!file.exists(gff_file)) {
        cat(sprintf("  MISSING: C output file %s not found.\n", gff_file))
        return()
    }
    
    # Read C GFF output (skipping header lines starting with #)
    c_data <- read.table(gff_file, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)
    
    # GFF columns: SeqID, Source, Type, Start, End, Score, Strand, Phase, Attributes
    c_starts <- c_data[[4]]
    c_ends   <- c_data[[5]]
    
    r_starts <- r_df$start
    r_ends   <- r_df$end
    
    # Check if counts match
    if (length(c_starts) == length(r_starts)) {
        cat(sprintf("  SUCCESS: Both found %d motifs.\n", length(r_starts)))
    } else {
        cat(sprintf("  FAILURE: Count mismatch! C=%d, R=%d\n", length(c_starts), length(r_starts)))
    }
    
    # Check if coordinates match
    matches <- sum(r_starts %in% c_starts & r_ends %in% c_ends)
    if (matches == length(c_starts) && length(c_starts) == length(r_starts)) {
        cat("  SUCCESS: All coordinates match perfectly.\n")
    } else {
        cat(sprintf("  WARNING: %d/%d coordinates match.\n", matches, max(length(c_starts), length(r_starts))))
    }
}

# 5. Execute comparisons
compare_results(r_results$apr, paste0(output_prefix, "_APR.gff"), "A-Phased Repeats (APR)")
compare_results(r_results$dr,  paste0(output_prefix, "_DR.gff"),  "Direct Repeats (DR)")
compare_results(r_results$gq,  paste0(output_prefix, "_GQ.gff"),  "G-quadruplexes (GQ)")
compare_results(r_results$ir,  paste0(output_prefix, "_IR.gff"),  "Inverted Repeats (IR)")
compare_results(r_results$mr,  paste0(output_prefix, "_MR.gff"),  "Mirror Repeats (MR)")
compare_results(r_results$str, paste0(output_prefix, "_STR.gff"), "Short Tandem Repeats (STR)")
compare_results(r_results$zdna, paste0(output_prefix, "_Z.gff"),   "Z-DNA Motifs")
