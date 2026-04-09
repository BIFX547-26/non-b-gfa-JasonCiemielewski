#' Find A-Phased Repeats (APR)
#'
#' @description
#' This function identifies A-Phased Repeats (APRs) in a DNA sequence. APRs are
#' associated with bent DNA structures and are characterized by runs of
#' Adenine (A) or Thymine (T) bases with a specific periodicity. The function
#' internally searches for both A-tracts and T-tracts by analyzing both the
#' forward and reverse-complement strands.
#'
#' @param dna_seq A character string representing the DNA sequence to analyze.
#' @param minAPR Minimum length of an A-tract (default: 3).
#' @param maxAPR Maximum length of an A-tract (default: 9).
#' @param minATracts Minimum number of consecutive A-tracts required to define an APR (default: 3).
#'
#' @return A data.frame containing the detected APRs with columns:
#' \itemize{
#'   \item \code{start}: 1-based start position of the motif.
#'   \item \code{loop}: Not used for APR (default: 0).
#'   \item \code{len}: Number of A-tracts in the repeat.
#'   \item \code{num}: Number of A-tracts in the repeat (redundant with len).
#'   \item \code{end}: 1-based end position of the motif.
#'   \item \code{sub}: Not used for APR (default: 0).
#'   \item \code{strand}: DNA strand (0 for plus).
#'   \item \code{special}: Internal classification flag.
#' }
#'
#' @examples
#' # Example 1: Sequence containing an APR
#' find_apr("ttttgaagaggaaatgaagggtattN")
#' 
#' # Example 2: Sequence with no APRs
#' find_apr("GCGCGCGCGCGCGCGCGCGC")
#' @export
find_apr <- function(dna_seq, minAPR = 3, maxAPR = 9, minATracts = 3) {
    # Input validation: Ensure sequence is character
    if (!is.character(dna_seq)) stop("dna_seq must be a character string.")
    
    # Call the Rcpp implementation for performance
    find_apr_rcpp(dna_seq, minAPR, maxAPR, minATracts)
}

#' Find Direct Repeats (DR)
#'
#' @description
#' This function identifies Direct Repeats (DRs) in a DNA sequence. DRs are 
#' associated with slipped-strand DNA structures and consist of two or more 
#' identical or nearly identical sequences in the same orientation.
#'
#' @param dna_seq A character string representing the DNA sequence to analyze.
#' @param mindir Minimum size of the repeating unit (default: 10).
#' @param maxdir Maximum size of the repeating unit (default: 300).
#' @param dspacer Maximum spacer length between repeating units (default: 10).
#'
#' @return A data.frame containing the detected DRs with columns:
#' \itemize{
#'   \item \code{start}: 1-based start position of the motif.
#'   \item \code{loop}: Spacer length between repeating units.
#'   \item \code{len}: Length of the repeating unit.
#'   \item \code{num}: Number of full repeating units found.
#'   \item \code{end}: 1-based end position of the motif.
#'   \item \code{sub}: Number of additional bases matching the repeating unit (partial repeat).
#'   \item \code{strand}: DNA strand (0 for plus).
#'   \item \code{special}: Internal classification flag (1 if classified as Slipped DNA).
#' }
#'
#' @examples
#' # Example 1: Sequence containing a Direct Repeat
#' find_dr("ACCCCTACCCCTACCCCTACCCCN", mindir = 5)
#' 
#' # Example 2: Sequence with no Direct Repeats
#' find_dr("ATGCATGCATGC", mindir = 10)
#' @export
find_dr <- function(dna_seq, mindir = 10, maxdir = 300, dspacer = 10) {
    # Input validation: Ensure sequence is character
    if (!is.character(dna_seq)) stop("dna_seq must be a character string.")
    
    # Call the Rcpp implementation for performance
    find_dr_rcpp(dna_seq, mindir, maxdir, dspacer)
}

#' Find G-quadruplexes (GQ)
#'
#' @description
#' This function identifies potential G-quadruplex (GQ) motifs in a DNA sequence.
#' G-quadruplexes are four-stranded DNA structures formed by G-rich sequences.
#'
#' @param dna_seq A character string representing the DNA sequence to analyze.
#' @param minGQ Minimum size of the Guanine (G) runs (default: 3).
#' @param maxGQspacer Maximum spacer (loop) length between G-runs (default: 7).
#'
#' @return A data.frame containing the detected GQs with columns:
#' \itemize{
#'   \item \code{start}: 1-based start position of the motif.
#'   \item \code{loop}: Not used for GQ (default: 0).
#'   \item \code{len}: Maximum G-run size found within the motif.
#'   \item \code{num}: Number of G-runs found in the motif.
#'   \item \code{end}: 1-based end position of the motif.
#'   \item \code{sub}: Number of consecutive G-islands identified.
#'   \item \code{strand}: DNA strand (0 for plus, 1 for minus).
#'   \item \code{special}: Internal classification flag.
#' }
#'
#' @examples
#' # Example 1: Sequence containing a G-quadruplex
#' find_gq("GGGAAGGAGGGAGAGACGGGAGGGN")
#' 
#' # Example 2: Sequence with no G-quadruplexes
#' find_gq("ATATATATATATATATATAT")
#' @export
find_gq <- function(dna_seq, minGQ = 3, maxGQspacer = 7) {
    # Input validation: Ensure sequence is character
    if (!is.character(dna_seq)) stop("dna_seq must be a character string.")
    
    # Call the Rcpp implementation for performance
    find_gq_rcpp(dna_seq, minGQ, maxGQspacer)
}

#' Find Inverted Repeats (IR)
#'
#' @description
#' This function identifies Inverted Repeats (IRs) in a DNA sequence. IRs are 
#' associated with cruciform DNA structures and consist of two sequences that 
#' are reverse complements of each other.
#'
#' @param dna_seq A character string representing the DNA sequence to analyze.
#' @param mincrf Minimum size of the repeating unit (default: 6).
#' @param cspacer Maximum spacer length between inverted units (default: 100).
#' @param cut Stem size threshold for short spacer constraint (default: 9).
#' @param shortSpacer Maximum spacer for short stems (below 'cut' threshold) (default: 4).
#'
#' @return A data.frame containing the detected IRs with columns:
#' \itemize{
#'   \item \code{start}: 1-based start position of the motif.
#'   \item \code{loop}: Spacer length between the inverted units.
#'   \item \code{len}: Length of the repeating unit (stem length).
#'   \item \code{num}: Number of permuations/occurrences (typically 1).
#'   \item \code{end}: 1-based end position of the motif.
#'   \item \code{sub}: Minimum loop boundary position.
#'   \item \code{strand}: DNA strand (0 for plus).
#'   \item \code{special}: Internal classification flag (1 if classified as Cruciform).
#' }
#'
#' @examples
#' # Example 1: Sequence containing an Inverted Repeat
#' find_ir("GCGCGCGCGCATATGCGCGCGCGC")
#' 
#' # Example 2: Sequence with no Inverted Repeats
#' find_ir("AAAAAAAAAAAAAAAAAAAA")
#' @export
find_ir <- function(dna_seq, mincrf = 6, cspacer = 100, cut = 9, shortSpacer = 4) {
    # Input validation: Ensure sequence is character
    if (!is.character(dna_seq)) stop("dna_seq must be a character string.")
    
    # Call the Rcpp implementation for performance
    find_ir_rcpp(dna_seq, mincrf, cspacer, cut, shortSpacer)
}

#' Find Mirror Repeats (MR)
#'
#' @description
#' This function identifies Mirror Repeats (MRs) in a DNA sequence. MRs are 
#' associated with triplex DNA structures and consist of two identical 
#' sequences in reverse orientation on the same strand.
#'
#' @param dna_seq A character string representing the DNA sequence to analyze.
#' @param minmir Minimum size of the repeating unit (default: 10).
#' @param mspacer Maximum spacer length between mirror units (default: 100).
#'
#' @return A data.frame containing the detected MRs with columns:
#' \itemize{
#'   \item \code{start}: 1-based start position of the motif.
#'   \item \code{loop}: Spacer length between mirror units.
#'   \item \code{len}: Length of the repeating unit.
#'   \item \code{num}: Number of permutations/occurrences.
#'   \item \code{end}: 1-based end position of the motif.
#'   \item \code{sub}: Minimum loop boundary position.
#'   \item \code{strand}: DNA strand (0 for plus).
#'   \item \code{special}: Internal classification flag (1 if classified as Triplex).
#' }
#'
#' @examples
#' # Example 1: Sequence containing a Mirror Repeat
#' find_mr("CCCCGCCGCCCGCCCGCCGCCCCN")
#' 
#' # Example 2: Sequence with no Mirror Repeats
#' find_mr("ATGCATGCATGCATGC")
#' @export
find_mr <- function(dna_seq, minmir = 10, mspacer = 100) {
    # Input validation: Ensure sequence is character
    if (!is.character(dna_seq)) stop("dna_seq must be a character string.")
    
    # Call the Rcpp implementation for performance
    find_mr_rcpp(dna_seq, minmir, mspacer)
}

#' Find Short Tandem Repeats (STR)
#'
#' @description
#' This function identifies Short Tandem Repeats (STRs) or microsatellites in 
#' a DNA sequence. STRs consist of short repeating units (1-9 bp) occurring 
#' multiple times in tandem.
#'
#' @param dna_seq A character string representing the DNA sequence to analyze.
#' @param minSTR Minimum unit size (default: 1).
#' @param maxSTR Maximum unit size (default: 9).
#' @param minSTRlen Minimum total length of the repeat region (default: 10).
#' @param minReps Minimum number of repeats required (default: 3).
#'
#' @return A data.frame containing the detected STRs with columns:
#' \itemize{
#'   \item \code{start}: 1-based start position of the motif.
#'   \item \code{loop}: Classification code for the STR type.
#'   \item \code{len}: Unit size of the repeat.
#'   \item \code{num}: Number of full repeating units.
#'   \item \code{end}: 1-based end position of the motif.
#'   \item \code{sub}: Number of additional partial matching bases at the end.
#'   \item \code{strand}: DNA strand (0 for plus).
#'   \item \code{special}: Internal classification flag.
#' }
#'
#' @examples
#' # Example 1: Sequence containing an STR (GT repeats)
#' find_str("GTGTGTGTGTGT")
#' 
#' # Example 2: Sequence with no STRs
#' find_str("ATGCATGCATGC", minReps = 10)
#' @export
find_str <- function(dna_seq, minSTR = 1, maxSTR = 9, minSTRlen = 10, minReps = 3) {
    # Input validation: Ensure sequence is character
    if (!is.character(dna_seq)) stop("dna_seq must be a character string.")
    
    # Call the Rcpp implementation for performance
    find_str_rcpp(dna_seq, minSTR, maxSTR, minSTRlen, minReps)
}

#' Find Z-DNA Motifs
#'
#' @description
#' This function identifies potential Z-DNA motifs in a DNA sequence. Z-DNA 
#' is a left-handed helical structure favored by alternating purine-pyrimidine 
#' sequences, especially (CG)n.
#'
#' @param dna_seq A character string representing the DNA sequence to analyze.
#' @param minZ Minimum run size required for detection (default: 10).
#'
#' @return A data.frame containing the detected Z-DNA motifs with columns:
#' \itemize{
#'   \item \code{start}: 1-based start position of the motif.
#'   \item \code{loop}: Score calculated based on purine-pyrimidine patterns.
#'   \item \code{len}: Length of the motif.
#'   \item \code{num}: Not used for Z-DNA (default: 0).
#'   \item \code{end}: 1-based end position of the motif.
#'   \item \code{sub}: Not used for Z-DNA (default: 0).
#'   \item \code{strand}: DNA strand (0 for plus).
#'   \item \code{special}: Internal classification flag (1 if score above threshold).
#' }
#'
#' @examples
#' # Example 1: Sequence containing a Z-DNA motif
#' find_zdna("gtgcacgcgtgcgtgN")
#' 
#' # Example 2: Sequence with no Z-DNA motifs
#' find_zdna("AAAAAAAAAAAAAA")
#' @export
find_zdna <- function(dna_seq, minZ = 10) {
    # Input validation: Ensure sequence is character
    if (!is.character(dna_seq)) stop("dna_seq must be a character string.")
    
    # Call the Rcpp implementation for performance
    find_zdna_rcpp(dna_seq, minZ)
}
