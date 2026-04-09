#' Analyze all non-B DNA motifs in a sequence
#' 
#' @description
#' This is a master function that runs all available motif-finding algorithms 
#' on a single DNA sequence. It consolidates results for APR, DR, GQ, IR, MR, 
#' STR, and Z-DNA into a structured list.
#'
#' @param dna_seq A character string representing the DNA sequence to analyze.
#' @return A named list of data.frames, where each element corresponds to a 
#' motif type:
#' \itemize{
#'   \item \code{apr}: A-Phased Repeats.
#'   \item \code{dr}: Direct Repeats (Slipped DNA).
#'   \item \code{gq}: G-quadruplexes.
#'   \item \code{ir}: Inverted Repeats (Cruciforms).
#'   \item \code{mr}: Mirror Repeats (Triplexes).
#'   \item \code{str}: Short Tandem Repeats.
#'   \item \code{zdna}: Z-DNA motifs.
#' }
#' 
#' @examples
#' seq <- "AGTGCAACCCAGAGGGCAGGATTTCCTGCTGGACTTTGAAATCCAACCCGGTCACCTACCCGCGCGACTG"
#' all_motifs <- gfa_analyze(seq)
#' names(all_motifs)
#' head(all_motifs$gq)
#' @export
gfa_analyze <- function(dna_seq) {
    # Input validation: Ensure sequence is character
    if (!is.character(dna_seq)) stop("dna_seq must be a character string.")

    # Execute all motif-finding functions and consolidate into a list
    list(
        apr = find_apr(dna_seq),
        dr  = find_dr(dna_seq),
        gq  = find_gq(dna_seq),
        ir  = find_ir(dna_seq),
        mr  = find_mr(dna_seq),
        str = find_str(dna_seq),
        zdna = find_zdna(dna_seq)
    )
}
