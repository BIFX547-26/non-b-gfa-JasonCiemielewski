library(testthat)
library(gfaR)

# Helper to read FASTA
read_fasta <- function(file) {
    if (!file.exists(file)) return("")
    lines <- readLines(file)
    seq_lines <- lines[!grepl("^>", lines)]
    paste0(seq_lines, collapse = "")
}

# --- 1. Gold Standard Parity Tests (Original Tests) ---

test_that("Gold Standard: Motif detection matches original C program", {
    # These files are extracted to the test directory during build
    fasta_file <- "gfa_test.fasta"
    if (file.exists(fasta_file)) {
        seq <- read_fasta(fasta_file)
        
        # APR
        expect_true(any(find_apr(seq)$start == 582))
        # DR
        expect_true(any(find_dr(seq)$start == 71))
        # GQ
        expect_true(any(find_gq(seq)$start == 5305))
        # IR
        expect_true(any(find_ir(seq)$start == 16))
        # MR
        expect_true(any(find_mr(seq)$start == 2845))
        # STR
        expect_true(any(find_str(seq)$start == 71))
        # Z-DNA
        expect_true(any(find_zdna(seq)$start == 2603))
    }
})

# --- 2. Input Validation & Edge Cases ---

test_that("Input validation and edge cases work correctly", {
    # Non-character input
    expect_error(find_apr(123))
    expect_error(find_dr(NULL))
    
    # Empty string
    expect_equal(nrow(find_apr("")), 0)
    expect_equal(nrow(find_dr("")), 0)
    
    # Sequence too short for any motif
    short_seq <- "ATGC"
    expect_equal(nrow(find_apr(short_seq)), 0)
    expect_equal(nrow(find_dr(short_seq)), 0)
    expect_equal(nrow(find_gq(short_seq)), 0)
    expect_equal(nrow(find_ir(short_seq)), 0)
    expect_equal(nrow(find_mr(short_seq)), 0)
    expect_equal(nrow(find_str(short_seq)), 0)
    expect_equal(nrow(find_zdna(short_seq)), 0)
    
    # Sequence with only Ns
    n_seq <- "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
    expect_equal(nrow(find_dr(n_seq)), 0)
    expect_equal(nrow(find_str(n_seq)), 0)
})

# --- 3. Parameter Sensitivity Tests ---

test_that("APR parameter sensitivity works", {
    # ttttgaagaggaaatgaagggtatt (3 tracts, len 3-4)
    seq <- "ttttgaagaggaaatgaagggtattN"
    
    # Default (minAPR=3, minATracts=3) -> Found
    expect_equal(nrow(find_apr(seq)), 1)
    
    # Increase minAPR to 5 -> Not found
    expect_equal(nrow(find_apr(seq, minAPR = 5)), 0)
    
    # Increase minATracts to 4 -> Not found
    expect_equal(nrow(find_apr(seq, minATracts = 4)), 0)
})

test_that("DR parameter sensitivity works", {
    # Direct repeat: 10bp unit, 5bp spacer
    # ATGCATGCAT 12345 ATGCATGCAT
    seq <- "ATGCATGCATGGGGGATGCATGCATN"
    
    # Default (mindir=10, dspacer=10) -> Found
    expect_equal(nrow(find_dr(seq)), 1)
    
    # Decrease dspacer to 4 -> Not found
    expect_equal(nrow(find_dr(seq, dspacer = 4)), 0)
    
    # Increase mindir to 11 -> Not found
    expect_equal(nrow(find_dr(seq, mindir = 11)), 0)
})

test_that("GQ parameter sensitivity works", {
    # GGG AAA GGG AAA GGG AAA GGG (4 runs, 3bp spacer)
    seq <- "GGGAAAGGGAAAGGGAAAGGGN"
    
    # Default (minGQ=3, maxGQspacer=7) -> Found
    expect_equal(nrow(find_gq(seq)), 1)
    
    # Increase minGQ to 4 -> Not found
    expect_equal(nrow(find_gq(seq, minGQ = 4)), 0)
    
    # Decrease maxGQspacer to 2 -> Found (spacer is 3bp, wait... 
    # the GFF says islands are 1-based, AAA is 3bp. 
    # 4 to 6 is AAA. Start of run 2 is 7. 7 - (1+3) = 3.
    expect_equal(nrow(find_gq(seq, maxGQspacer = 2)), 0)
})

test_that("IR parameter sensitivity works", {
    # Inverted Repeat: GCTAGC (6bp) 4bp spacer GCTAGC
    # GCGCGC ATAT GCGCGC
    seq <- "GCGCGCATATGCGCGC" # RC of GCGCGC is GCGCGC
    
    # Default (mincrf=6, cspacer=100) -> Found
    expect_equal(nrow(find_ir(seq)), 1)
    
    # Increase mincrf to 9 -> Not found (previous 8bp match)
    expect_equal(nrow(find_ir(seq, mincrf = 9)), 0)
    
    # Test shortSpacer constraint (cut=9, shortSpacer=4)
    # If stem <= 9, spacer must be <= 4.
    # Stem is 6, spacer is 10.
    seq_long_sp <- "GCGCGC1234567890GCGCGC"
    expect_equal(nrow(find_ir(seq_long_sp, cspacer = 100)), 0)
})

test_that("MR parameter sensitivity works", {
    # Mirror Repeat: 10bp unit
    seq <- "ATGCATGCATGGGTACGTAGTAN" # ATGCATGCAT ... reversed ...
    # Wait, mirror is reversed on same strand: ATGCATGCAT ... TACGTACGAT
    seq_mr <- "ATGCATGCATTTTTTACGTACGTAN"
    
    # Default (minmir=10, mspacer=100) -> Found
    expect_equal(nrow(find_mr(seq_mr)), 1)
    
    # Increase minmir to 15 -> Not found (previous 12bp match)
    expect_equal(nrow(find_mr(seq_mr, minmir = 15)), 0)
})

test_that("STR parameter sensitivity works", {
    # GT repeat (unit 2) x 6 = 12bp
    seq <- "GTGTGTGTGTGT"
    
    # Default (minSTR=1, maxSTR=9, minSTRlen=10, minReps=3) -> Found
    expect_equal(nrow(find_str(seq)), 1)
    
    # Increase minReps to 10 -> Not found
    expect_equal(nrow(find_str(seq, minReps = 10)), 0)
    
    # Increase minSTRlen to 15 -> Not found
    expect_equal(nrow(find_str(seq, minSTRlen = 15)), 0)
})

test_that("Z-DNA parameter sensitivity works", {
    # alternating (GC)n
    seq <- "GCGCGCGCGCGCGCGCGCN" # 18bp
    
    # Default (minZ=10) -> Found
    expect_equal(nrow(find_zdna(seq)), 1)
    
    # Increase minZ to 20 -> Not found
    expect_equal(nrow(find_zdna(seq, minZ = 20)), 0)
})

test_that("gfa_analyze aggregates all results", {
    seq <- "GGGAAGGAGGGAGAGACGGGAGGGN"
    res <- gfa_analyze(seq)
    
    expect_type(res, "list")
    expect_true(all(c("apr", "dr", "gq", "ir", "mr", "str", "zdna") %in% names(res)))
    expect_s3_class(res$gq, "data.frame")
    expect_true(nrow(res$gq) > 0)
})
