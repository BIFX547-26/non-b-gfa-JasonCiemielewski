# Comparison script: GFA C-program vs gfaR R package
library(gfaR)

# 1. Setup paths
fasta_file <- "test_data/gfa_test.fasta"
c_src_dir  <- "gfa_c_program"
c_program  <- file.path(c_src_dir, "gfa.exe")
output_prefix <- "test_data/c_output"

# 2. Auto-compile C program if missing
if (!file.exists(c_program)) {
    cat("--- Compiling original C program ---\n")
    c_files <- c("gfa.c", "cdna.c", "findAPR.c", "findDR.c", "findGQ.c", 
                 "findIR.c", "findMR.c", "findSTR.c", "findZDNA.c", 
                 "is_subset.c", "nulls.c", "print_gff_file.c", 
                 "print_tsv_file.c", "print_usage.c", "process_repeats.c", 
                 "rcdna.c", "read_fasta.c", "read_mult_fasta.c")
    
    # Attempt to find gcc (common in Rtools on Windows)
    gcc_path <- "gcc"
    # Check if we are on Windows and need to look for Rtools
    if (.Platform$OS.type == "windows") {
        # This is a common path for Rtools45, matching your environment
        gcc_hint <- "C:/rtools45/x86_64-w64-mingw32.static.posix/bin/gcc.exe"
        if (file.exists(gcc_hint)) gcc_path <- gcc_hint
    }
    
    compile_cmd <- sprintf("%s -O2 %s -o %s -lm", 
                           gcc_path, 
                           paste(file.path(c_src_dir, c_files), collapse = " "), 
                           c_program)
    
    system(compile_cmd)
    
    if (!file.exists(c_program)) {
        stop("Failed to compile original C program. Please ensure a C compiler (gcc) is available.")
    }
}

# 3. Run the original C program
cat("--- Running original C program ---\n")
system(paste(c_program, "-seq", fasta_file, "-chrom Test -out", output_prefix))

# 4. Run the R package
cat("\n--- Running gfaR R package ---\n")
read_fasta <- function(file) {
    lines <- readLines(file)
    seq_lines <- lines[!grepl("^>", lines)]
    paste0(seq_lines, collapse = "")
}
dna_seq <- read_fasta(fasta_file)
r_results <- gfa_analyze(dna_seq)

# 5. Define enhanced comparison function
compare_results <- function(r_df, gff_file, motif_name) {
    cat(sprintf("\n========================================================\n"))
    cat(sprintf(" MOTIF TYPE: %s\n", motif_name))
    cat(sprintf("========================================================\n"))
    
    if (!file.exists(gff_file)) {
        cat(sprintf("  MISSING: C output file %s not found.\n", gff_file))
        return()
    }
    
    # Read C GFF output
    c_data <- read.table(gff_file, comment.char = "#", sep = "\t", stringsAsFactors = FALSE)
    c_comp <- data.frame(Version = "C", start = c_data[[4]], end = c_data[[5]])
    
    # Prepare R data for comparison
    r_comp <- data.frame(Version = "R", start = r_df$start, end = r_df$end)
    
    # Show counts
    cat(sprintf(" Counts Found -> C: %d | R: %d\n\n", nrow(c_comp), nrow(r_comp)))
    
    # Create side-by-side view
    max_rows <- min(10, nrow(c_comp), nrow(r_comp))
    
    if (max_rows > 0) {
        cat(" Side-by-Side Comparison (First 10 matches):\n")
        comparison_table <- data.frame(
            Index = 1:max_rows,
            C_Start = c_comp$start[1:max_rows],
            R_Start = r_comp$start[1:max_rows],
            C_End   = c_comp$end[1:max_rows],
            R_End   = r_comp$end[1:max_rows]
        )
        comparison_table$Match <- ifelse(
            comparison_table$C_Start == comparison_table$R_Start & 
            comparison_table$C_End == comparison_table$R_End, 
            "YES", "NO"
        )
        print(comparison_table, row.names = FALSE)
    }
    
    # Summary Statement
    matches <- sum(r_df$start %in% c_comp$start & r_df$end %in% c_comp$end)
    if (matches == nrow(c_comp) && nrow(c_comp) == nrow(r_comp)) {
        cat("\n RESULT: SUCCESS - All records are identical.\n")
    } else {
        cat(sprintf("\n RESULT: DISCREPANCY - Found %d coordinate matches.\n", matches))
    }
}

# 6. Execute comparisons
compare_results(r_results$apr, paste0(output_prefix, "_APR.gff"), "A-Phased Repeats (APR)")
compare_results(r_results$dr,  paste0(output_prefix, "_DR.gff"),  "Direct Repeats (DR)")
compare_results(r_results$gq,  paste0(output_prefix, "_GQ.gff"),  "G-quadruplexes (GQ)")
compare_results(r_results$ir,  paste0(output_prefix, "_IR.gff"),  "Inverted Repeats (IR)")
compare_results(r_results$mr,  paste0(output_prefix, "_MR.gff"),  "Mirror Repeats (MR)")
compare_results(r_results$str, paste0(output_prefix, "_STR.gff"), "Short Tandem Repeats (STR)")
compare_results(r_results$zdna, paste0(output_prefix, "_Z.gff"),   "Z-DNA Motifs")
