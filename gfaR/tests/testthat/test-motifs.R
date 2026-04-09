library(testthat)
library(gfaR)

read_fasta <- function(file) {
    lines <- readLines(file)
    seq_lines <- lines[!grepl("^>", lines)]
    paste0(seq_lines, collapse = "")
}

test_that("APR detection matches original C program", {
    fasta_file <- "gfa_test.fasta"
    seq <- read_fasta(fasta_file)
    res <- find_apr(seq)
    expect_true(any(res$start == 582 & res$end == 606 & res$num == 3))
})

test_that("DR detection matches original C program", {
    fasta_file <- "gfa_test.fasta"
    seq <- read_fasta(fasta_file)
    res <- find_dr(seq)
    # Test 71 93
    expect_true(any(res$start == 71 & res$end == 93 & res$loop == 1))
})

test_that("GQ detection matches original C program", {
    fasta_file <- "gfa_test.fasta"
    seq <- read_fasta(fasta_file)
    res <- find_gq(seq)
    # Test 5305 5328
    expect_true(any(res$start == 5305 & res$end == 5328))
})

test_that("IR detection matches original C program", {
    fasta_file <- "gfa_test.fasta"
    seq <- read_fasta(fasta_file)
    res <- find_ir(seq)
    # Test 16 29
    expect_true(any(res$start == 16 & res$end == 29))
})

test_that("MR detection matches original C program", {
    fasta_file <- "gfa_test.fasta"
    seq <- read_fasta(fasta_file)
    res <- find_mr(seq)
    # Test 2845 2867
    expect_true(any(res$start == 2845 & res$end == 2867))
})

test_that("STR detection matches original C program", {
    fasta_file <- "gfa_test.fasta"
    seq <- read_fasta(fasta_file)
    res <- find_str(seq)
    # Test 71 93
    expect_true(any(res$start == 71 & res$end == 93))
})

test_that("Z-DNA detection matches original C program", {
    fasta_file <- "gfa_test.fasta"
    seq <- read_fasta(fasta_file)
    res <- find_zdna(seq)
    # Test 2603 2614
    expect_true(any(res$start == 2603 & res$end == 2614 & res$loop == 16))
})
