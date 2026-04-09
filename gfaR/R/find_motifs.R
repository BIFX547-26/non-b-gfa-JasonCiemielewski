#' Find A-Phased Repeats
#' @param dna_seq DNA sequence string
#' @param minAPR Minimum length of A-tract (default: 3)
#' @param maxAPR Maximum length of A-tract (default: 9)
#' @param minATracts Minimum number of A-tracts (default: 3)
#' @export
find_apr <- function(dna_seq, minAPR = 3, maxAPR = 9, minATracts = 3) {
    find_apr_rcpp(dna_seq, minAPR, maxAPR, minATracts)
}

#' Find Direct Repeats
#' @param dna_seq DNA sequence string
#' @param mindir Minimum repeat size (default: 10)
#' @param maxdir Maximum repeat size (default: 300)
#' @param dspacer Maximum spacer size (default: 10)
#' @export
find_dr <- function(dna_seq, mindir = 10, maxdir = 300, dspacer = 10) {
    find_dr_rcpp(dna_seq, mindir, maxdir, dspacer)
}

#' Find G-quadruplexes
#' @param dna_seq DNA sequence string
#' @param minGQ Minimum G-run size (default: 3)
#' @param maxGQspacer Maximum loop size (default: 7)
#' @export
find_gq <- function(dna_seq, minGQ = 3, maxGQspacer = 7) {
    find_gq_rcpp(dna_seq, minGQ, maxGQspacer)
}

#' Find Inverted Repeats
#' @param dna_seq DNA sequence string
#' @param mincrf Minimum repeat size (default: 6)
#' @param cspacer Maximum spacer size (default: 100)
#' @param cut Stem size threshold for short spacer check (default: 9)
#' @param shortSpacer Maximum spacer for short stems (default: 4)
#' @export
find_ir <- function(dna_seq, mincrf = 6, cspacer = 100, cut = 9, shortSpacer = 4) {
    find_ir_rcpp(dna_seq, mincrf, cspacer, cut, shortSpacer)
}

#' Find Mirror Repeats
#' @param dna_seq DNA sequence string
#' @param minmir Minimum repeat size (default: 10)
#' @param mspacer Maximum spacer size (default: 100)
#' @export
find_mr <- function(dna_seq, minmir = 10, mspacer = 100) {
    find_mr_rcpp(dna_seq, minmir, mspacer)
}

#' Find Short Tandem Repeats
#' @param dna_seq DNA sequence string
#' @param minSTR Minimum unit size (default: 1)
#' @param maxSTR Maximum unit size (default: 9)
#' @param minSTRlen Minimum total length (default: 10)
#' @param minReps Minimum number of repeats (default: 3)
#' @export
find_str <- function(dna_seq, minSTR = 1, maxSTR = 9, minSTRlen = 10, minReps = 3) {
    find_str_rcpp(dna_seq, minSTR, maxSTR, minSTRlen, minReps)
}

#' Find Z-DNA
#' @param dna_seq DNA sequence string
#' @param minZ Minimum run size (default: 10)
#' @export
find_zdna <- function(dna_seq, minZ = 10) {
    find_zdna_rcpp(dna_seq, minZ)
}
