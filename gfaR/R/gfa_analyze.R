#' Analyze all non-B DNA motifs in a sequence
#' 
#' @param dna_seq DNA sequence string
#' @return A named list of data.frames for each motif type
#' @export
gfa_analyze <- function(dna_seq) {
    list(
        apr = find_apr(dna_seq),
        dr = find_dr(dna_seq),
        gq = find_gq(dna_seq),
        ir = find_ir(dna_seq),
        mr = find_mr(dna_seq),
        str = find_str(dna_seq),
        zdna = find_zdna(dna_seq)
    )
}
