/*
 * find_motifs.cpp
 * 
 * This code is an Rcpp refactor of the original Genomic Feature Analyzer (GFA)
 * C-program developed by Cer et al. (2013).
 * 
 * Original Authors: Regina Z. Cer, Duncan E. Donohue, Uma S. Mudunuri, et al.
 * Refactored by: Jason Ciemielewski (2026)
 * 
 * Original Publication:
 * Cer, R. Z., et al. Nucl. Acids Res. (2013) 41 (D1): D94-D100.
 * doi: 10.1093/nar/gks955
 */

#include <Rcpp.h>
#include "gfa.h"
#include <iostream>
#include <cmath>

using namespace Rcpp;

/* Helper to convert REP vector to DataFrame */
DataFrame reps_to_df(const std::vector<REP>& reps) {
    int n = reps.size();
    IntegerVector start(n), loop(n), len(n), num(n), end(n), sub(n), strand(n), special(n);
    for (int i = 0; i < n; ++i) {
        start[i] = reps[i].start;
        loop[i] = reps[i].loop;
        len[i] = reps[i].len;
        num[i] = reps[i].num;
        end[i] = reps[i].end;
        sub[i] = reps[i].sub;
        strand[i] = reps[i].strand;
        special[i] = reps[i].special;
    }
    return DataFrame::create(
        _["start"] = start,
        _["loop"] = loop,
        _["len"] = len,
        _["num"] = num,
        _["end"] = end,
        _["sub"] = sub,
        _["strand"] = strand,
        _["special"] = special
    );
}

/* APR (A-Phased Repeats) Finder */
// [[Rcpp::export]]
DataFrame find_apr_rcpp(std::string dna_seq, int minAPR = 4, int maxAPR = 9, int minATracts = 3) {
    int total_bases = dna_seq.length();
    std::string dna2 = reverse_complement(dna_seq);
    
    std::vector<potential_Bent_DNA> pAPRs;
    int nAs = 0;
    int i = 0;
    
    while (i < total_bases) {
        char base = std::tolower(dna_seq[i]);
        if (base == 'a' || base == 't') {
            nAs++;
        } else {
            if (nAs >= minAPR && nAs <= maxAPR) {
                // Rcout << "Found potential A-tract candidate: len=" << nAs << " at " << i - nAs + 1 << std::endl;
                int strt = i - nAs + 1;
                int ATend = strt + nAs; 
                int maxATlen = 0, maxTlen = 0, Alen = 0, Tlen = 0, ATlen = 0, TAlen = 0, maxATend = 0;
                int maxATlen_rc = 0, maxTlen_rc = 0, Alen_rc = 0, Tlen_rc = 0, ATlen_rc = 0, TAlen_rc = 0, maxATend_rc = 0;
                
                int n_rc = total_bases - ATend;
                for (int n = strt - 1; n < ATend - 1; n++) {
                    n_rc++;
                    char dn = std::tolower(dna_seq[n]);
                    char d2n = std::tolower(dna2[n_rc]);
                    
                    if (dn == 'a') {
                        Tlen = 0; TAlen = 0;
                        if (n > 0 && std::tolower(dna_seq[n-1]) == 't') { Alen = 0; ATlen = 0; }
                        else { Alen++; ATlen++; }
                    }
                    if (d2n == 'a') {
                        Tlen_rc = 0; TAlen_rc = 0;
                        if (n_rc > 0 && std::tolower(dna2[n_rc-1]) == 't') { Alen_rc = 0; ATlen_rc = 0; }
                        else { Alen_rc++; ATlen_rc++; }
                    }
                    if (dn == 't') {
                        if (TAlen < Alen) { TAlen++; ATlen++; }
                        else { Tlen++; TAlen = 0; ATlen = 0; Alen = 0; }
                    }
                    if (d2n == 't') {
                        if (TAlen_rc < Alen_rc) { TAlen_rc++; ATlen_rc++; }
                        else { Tlen_rc++; TAlen_rc = 0; ATlen_rc = 0; Alen_rc = 0; }
                    }
                    if (maxATlen < ATlen) { maxATlen = ATlen; maxATend = n; }
                    if (maxTlen < Tlen) { maxTlen = Tlen; }
                    if (maxATlen_rc < ATlen_rc) { maxATlen_rc = ATlen_rc; maxATend_rc = n_rc; }
                    if (maxTlen_rc < Tlen_rc) { maxTlen_rc = Tlen_rc; }
                }
                
                if (((maxATlen - maxTlen) >= minAPR) || ((maxATlen_rc - maxTlen_rc) >= minAPR)) {
                    potential_Bent_DNA p;
                    p.end = strt + nAs;
                    p.strt = strt;
                    if ((maxATlen - maxTlen) >= (maxATlen_rc - maxATlen_rc)) {
                        p.a_center = ((double) maxATend - (((double) maxATlen - 1) / 2)) + 1;
                    } else {
                        p.a_center = total_bases - (((double) maxATend_rc - (((double) maxATlen_rc - 1) / 2)));
                    }
                    pAPRs.push_back(p);
                }

            }
            nAs = 0;
        }
        i++;
    }
    // Rcout << "nProcessedATs: " << pAPRs.size() << std::endl;

    std::vector<REP> arep;
    int tracts = 1;
    double distToNext = 0;
    int nProcessedATs = pAPRs.size();
    
    for (int i = 0; i < nProcessedATs - 1; i++) {
        distToNext = pAPRs[i+1].a_center - pAPRs[i].a_center;
        // if (pAPRs[i].strt > 570 && pAPRs[i].strt < 620) Rcout << "distToNext: " << distToNext << " at " << pAPRs[i].strt << std::endl;
        if (distToNext <= 11.1 && distToNext >= 9.9) {
            tracts++;
        } else {
            if (tracts >= minATracts) {
                REP r;
                r.start = pAPRs[(i - tracts) + 1].strt;
                r.loop = 0;
                r.num = tracts;
                r.strand = 0;
                r.len = tracts;
                r.end = pAPRs[i].end - 1;
                r.sub = 0; r.special = 0;
                arep.push_back(r);
            }
            tracts = 1;
        }
    }
    // Handle final tracts
    if (tracts >= minATracts && nProcessedATs >= minATracts) {
        REP r;
        r.start = pAPRs[nProcessedATs - tracts].strt;
        r.loop = 0;
        r.num = tracts;
        r.strand = 0;
        r.len = tracts;
        r.end = pAPRs[nProcessedATs - 1].end - 1;
        r.sub = 0; r.special = 0;
        arep.push_back(r);
    }

    return reps_to_df(arep);
}

/* Direct Repeat Finder */
// [[Rcpp::export]]
DataFrame find_dr_rcpp(std::string dna_seq, int mindir = 10, int maxdir = 100, int dspacer = 100) {
    int total_bases = dna_seq.length();
    std::vector<REP> drep;
    int end = 0;
    int lasti = total_bases - (mindir * 2);

    for (int strti = 0; strti <= lasti; strti++) {
        while (strti < total_bases && std::tolower(dna_seq[strti]) == 'n') strti++;
        if (strti >= lasti) break;

        for (int size = maxdir; size >= mindir; size--) {
            if (((size * 2) + dspacer) <= (end - strti)) {
                size = mindir - 1; // equivalent to breaking inner loop
                continue;
            }
            int spMin = std::max(0, ((end - strti) - (size * 2)) + 2);
            int spMax = std::min(dspacer, lasti - strti);
            for (int sp = spMin; sp <= spMax; sp++) {
                int j = strti + size + sp;
                int i = strti;
                int k = 0;
                while (k < size && j < total_bases && std::tolower(dna_seq[i]) == std::tolower(dna_seq[j]) && std::tolower(dna_seq[i]) != 'n') {
                    k++; j++; i++;
                }
                if (k == size) {
                    int totlen = k;
                    if (sp == 0) {
                        while (j < total_bases && std::tolower(dna_seq[i]) == std::tolower(dna_seq[j])) {
                            totlen++; j++; i++;
                        }
                    }
                    REP r;
                    r.start = strti + 1;
                    r.len = size;
                    r.loop = sp;
                    r.num = totlen / r.len;
                    r.end = j;
                    r.sub = (totlen % r.len);
                    r.strand = 0;
                    r.special = 0;
                    drep.push_back(r);
                    
                    end = j - 1;
                    sp = dspacer + 1;
                    size = mindir - 1;
                }
            }
        }
    }
    return reps_to_df(drep);
}

/* G-quadruplex Finder */
// [[Rcpp::export]]
DataFrame find_gq_rcpp(std::string dna_seq, int minGQ = 3, int maxGQspacer = 7) {
    int total_bases = dna_seq.length();
    std::vector<G_Island> gisle, rcgisle;
    int ngs = 0, ncs = 0;
    
    for (int i = 0; i <= total_bases; i++) {
        char base = (i < total_bases) ? std::tolower(dna_seq[i]) : '\0';
        if (base == 'g') {
            ngs++;
        } else {
            if (ngs >= minGQ) {
                G_Island g; g.strt = i - ngs + 1; g.len = ngs;
                gisle.push_back(g);
            }
            ngs = 0;
        }
        if (base == 'c') {
            ncs++;
        } else {
            if (ncs >= minGQ) {
                G_Island g; g.strt = i - ncs + 1; g.len = ncs;
                rcgisle.push_back(g);
            }
            ncs = 0;
        }
    }

    std::vector<REP> grep;
    for (int strand = 0; strand < 2; strand++) {
        std::vector<G_Island>& islands = (strand == 0) ? gisle : rcgisle;
        int nIls = islands.size();
        for (int i = 0; i < nIls; i++) {
            int conIls = 1;
            int npos = (int)(std::floor((islands[i].len + 1) / (minGQ + 1)));
            int i2 = i + 1;
            while (i2 < nIls && (islands[i2].strt - (islands[i2-1].strt + islands[i2-1].len)) <= maxGQspacer) {
                conIls++;
                npos += (int)(std::floor((islands[i2].len + 1) / (minGQ + 1)));
                i2++;
            }
            if (npos >= 4) {
                int maxGQ_found = minGQ;
                for (int j = i; j < i2; j++) {
                    for (int k = islands[j].len; k > maxGQ_found; k--) {
                        int nposMax = (int)(std::floor((islands[j].len + 1) / (k + 1)));
                        for (int m = j + 1; m < i2; m++) {
                            nposMax += (int)(std::floor((islands[m].len + 1) / (k + 1)));
                            if (nposMax >= 4) {
                                maxGQ_found = k;
                                break;
                            }
                        }
                        if (nposMax >= 4) break;
                    }
                }
                REP r;
                r.start = islands[i].strt;
                r.num = npos;
                r.sub = conIls;
                r.len = maxGQ_found;
                r.end = islands[i2-1].strt + islands[i2-1].len - 1;
                r.strand = strand;
                r.special = 0;
                r.loop = 0; // Default
                grep.push_back(r);
            }
            i = i + conIls - 1;
        }
    }
    return reps_to_df(grep);
}

/* Inverted Repeat Finder */
// [[Rcpp::export]]
DataFrame find_ir_rcpp(std::string dna_seq, int mincrf = 6, int cspacer = 100, int cut = 9, int shortSpacer = 4) {
    int total_bases = dna_seq.length();
    std::string dna3 = complement(dna_seq);
    std::vector<REP> irep;
    bool rightShifted = false;
    bool leftShifted = false;
    int maxcBack = 5;

    for (int strti = mincrf; strti <= (total_bases - mincrf); strti++) {
        int maxSP = std::min(cspacer, (total_bases - (strti + mincrf)));
        for (int sp = 0; sp <= maxSP; sp++) {
            int i = strti;
            int k = 0;
            int j = strti + sp + 1;
            while (j < total_bases && i >= 0 && std::tolower(dna_seq[i]) == std::tolower(dna3[j]) && std::tolower(dna_seq[j]) != 'n') {
                k++; j++; i--;
            }
            if (k >= mincrf) {
                if (k <= cut && sp > shortSpacer) continue;
                int tmpStart = (strti - k) + 2;
                int tmpStop = strti + k + sp + 1;
                
                if (irep.empty()) {
                    rightShifted = false; leftShifted = false;
                    REP r; r.start = tmpStart; r.sub = tmpStop; r.len = k; r.loop = sp; r.num = 1; r.end = tmpStop; r.strand = 0; r.special = 0;
                    irep.push_back(r);
                } else {
                    while (!irep.empty() && irep.back().end <= tmpStop && irep.back().start >= tmpStart && irep.back().len < k) {
                        irep.pop_back();
                        rightShifted = false; leftShifted = false;
                    }
                    while (!irep.empty() && irep.back().end >= tmpStop && irep.back().start <= tmpStart && irep.back().len < k) {
                        irep.pop_back();
                        rightShifted = false; leftShifted = false;
                    }
                    
                    if (!irep.empty() && irep.back().end <= tmpStop && irep.back().start >= tmpStart && irep.back().len > k) {
                        rightShifted = false; leftShifted = false;
                    } else if (!irep.empty() && irep.back().end >= tmpStop && irep.back().start <= tmpStart && irep.back().len > k) {
                        rightShifted = false; leftShifted = false;
                    } else if (!irep.empty() && tmpStop == irep.back().end && k == irep.back().len && !rightShifted) {
                        leftShifted = true; rightShifted = false;
                        irep.back().num++;
                        irep.back().sub = tmpStart;
                    } else if (!irep.empty() && tmpStart == irep.back().start && k == irep.back().len && !leftShifted) {
                        rightShifted = true; leftShifted = false;
                        if (irep.size() >= 2 && irep[irep.size()-2].end <= tmpStop && irep[irep.size()-2].start >= tmpStart && irep[irep.size()-2].len < k) {
                            irep[irep.size()-2] = irep.back();
                            irep.pop_back();
                        }
                        irep.back().num++;
                        irep.back().end = tmpStop;
                    } else {
                        rightShifted = false; leftShifted = false;
                        REP r; r.start = tmpStart; r.sub = tmpStop; r.len = k; r.loop = sp; r.num = 1; r.end = tmpStop; r.strand = 0; r.special = 0;
                        irep.push_back(r);
                        
                        for (int cBack = 1; cBack <= maxcBack; cBack++) {
                            int idx = irep.size() - 1;
                            int prevIdx = idx - cBack;
                            if (prevIdx < 0) break;
                            while (prevIdx >= 0 && (
                                (irep[prevIdx].end >= irep[idx].end && irep[prevIdx].start <= irep[idx].start) ||
                                (irep[prevIdx].end <= irep[idx].end && irep[prevIdx].start >= irep[idx].start)
                            )) {
                                if (irep[prevIdx].len == irep[idx].len) {
                                    if (irep[prevIdx].loop > irep[idx].loop) {
                                        irep.erase(irep.begin() + prevIdx);
                                        idx--;
                                    }
                                } else if (irep[prevIdx].len < irep[idx].len) {
                                    irep.erase(irep.begin() + prevIdx);
                                    idx--;
                                }
                                rightShifted = false; leftShifted = false;
                                irep.pop_back();
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    return reps_to_df(irep);
}

/* Mirror Repeat Finder */
// [[Rcpp::export]]
DataFrame find_mr_rcpp(std::string dna_seq, int minmir = 10, int mspacer = 100) {
    int total_bases = dna_seq.length();
    std::vector<REP> mrep;
    bool rightShifted = false;
    bool leftShifted = false;
    int maxcBack = 5;

    for (int strti = minmir; strti <= (total_bases - minmir); strti++) {
        int maxSP = std::min(mspacer, (total_bases - (strti + minmir)));
        for (int sp = 0; sp <= maxSP; sp++) {
            int i = strti;
            int k = 0;
            int j = strti + sp + 1;
            while (j < total_bases && i >= 0 && std::tolower(dna_seq[i]) == std::tolower(dna_seq[j]) && std::tolower(dna_seq[j]) != 'n') {
                k++; j++; i--;
            }
            if (k >= minmir) {
                int tmpStart = (strti - k) + 2;
                int tmpStop = strti + k + sp + 1;
                
                if (mrep.empty()) {
                    rightShifted = false; leftShifted = false;
                    REP r; r.start = tmpStart; r.sub = tmpStop; r.len = k; r.loop = sp; r.num = 1; r.end = tmpStop; r.strand = 0; r.special = 0;
                    mrep.push_back(r);
                } else {
                    // check for immediate inclusions
                    while (!mrep.empty() && mrep.back().end <= tmpStop && mrep.back().start >= tmpStart && mrep.back().len < k) {
                        mrep.pop_back();
                        rightShifted = false; leftShifted = false;
                    }
                    while (!mrep.empty() && mrep.back().end >= tmpStop && mrep.back().start <= tmpStart && mrep.back().len < k) {
                        mrep.pop_back();
                        rightShifted = false; leftShifted = false;
                    }
                    
                    if (!mrep.empty() && mrep.back().end <= tmpStop && mrep.back().start >= tmpStart && mrep.back().len > k) {
                        rightShifted = false; leftShifted = false;
                    } else if (!mrep.empty() && mrep.back().end >= tmpStop && mrep.back().start <= tmpStart && mrep.back().len > k) {
                        rightShifted = false; leftShifted = false;
                    } else if (!mrep.empty() && tmpStop == mrep.back().end && k == mrep.back().len && !rightShifted) {
                        leftShifted = true; rightShifted = false;
                        mrep.back().num++;
                        mrep.back().sub = tmpStart;
                    } else if (!mrep.empty() && tmpStart == mrep.back().start && k == mrep.back().len && !leftShifted) {
                        rightShifted = true; leftShifted = false;
                        if (mrep.size() >= 2 && mrep[mrep.size()-2].end <= tmpStop && mrep[mrep.size()-2].start >= tmpStart && mrep[mrep.size()-2].len < k) {
                            mrep[mrep.size()-2] = mrep.back();
                            mrep.pop_back();
                        }
                        mrep.back().num++;
                        mrep.back().end = tmpStop;
                    } else {
                        rightShifted = false; leftShifted = false;
                        REP r; r.start = tmpStart; r.sub = tmpStop; r.len = k; r.loop = sp; r.num = 1; r.end = tmpStop; r.strand = 0; r.special = 0;
                        mrep.push_back(r);
                        
                        // Back-check for overlaps
                        for (int cBack = 1; cBack <= maxcBack; cBack++) {
                            int idx = mrep.size() - 1;
                            int prevIdx = idx - cBack;
                            if (prevIdx < 0) break;
                            
                            while (prevIdx >= 0 && (
                                (mrep[prevIdx].end >= mrep[idx].end && mrep[prevIdx].start <= mrep[idx].start) ||
                                (mrep[prevIdx].end <= mrep[idx].end && mrep[prevIdx].start >= mrep[idx].start)
                            )) {
                                if (mrep[prevIdx].len == mrep[idx].len) {
                                    if (mrep[prevIdx].loop > mrep[idx].loop) {
                                        mrep.erase(mrep.begin() + prevIdx);
                                        idx--;
                                    } else {
                                        // current is worse? wait, C code deletes previous if its loop is larger
                                    }
                                } else if (mrep[prevIdx].len < mrep[idx].len) {
                                    mrep.erase(mrep.begin() + prevIdx);
                                    idx--;
                                }
                                rightShifted = false; leftShifted = false;
                                mrep.pop_back(); // Remove the current one because we erased previous? 
                                // WAIT! The C code --ndx; removes the CURRENT one!
                                break; 
                            }
                        }
                    }
                }
            }
        }
    }
    return reps_to_df(mrep);
}

/* Helper for STR classification */
int nonBstr(const std::string& dna, int start, int len) {
    int code = 0;
    bool isEven = (len % 2 == 0);
    bool isSymetric = true;
    bool isPUPY = true;
    bool isComp = true;

    if (len >= 2) {
        int j = start + len - 2;
        for (int i = 0; i <= (len / 2) - 1; i++) {
            if (std::tolower(dna[start + i - 1]) != std::tolower(dna[j])) isSymetric = false;
            if (isEven) {
                char c1 = std::tolower(dna[start + i - 1]);
                char c2 = std::tolower(dna[j]);
                if (c1 == 'a' && c2 != 't') isComp = false;
                else if (c1 == 't' && c2 != 'a') isComp = false;
                else if (c1 == 'c' && c2 != 'g') isComp = false;
                else if (c1 == 'g' && c2 != 'c') isComp = false;
            } else isComp = false;
            j--;
        }
        for (int i = start; i < (len + start - 1); i++) {
            char c1 = std::tolower(dna[i]);
            char c0 = std::tolower(dna[i-1]);
            bool r1 = (c1 == 'a' || c1 == 'g');
            bool r0 = (c0 == 'a' || c0 == 'g');
            bool y1 = (c1 == 't' || c1 == 'c');
            bool y0 = (c0 == 't' || c0 == 'c');
            if (r1 && r0) isPUPY = false;
            if (y1 && y0) isPUPY = false;
        }
    } else {
        isComp = false; isSymetric = true; isPUPY = false;
    }

    if (isEven) code += 1;
    if (isPUPY) code += 2;
    if (isSymetric) code += 4;
    if (isComp) code += 8;
    return code;
}

/* Short Tandem Repeat Finder */
// [[Rcpp::export]]
DataFrame find_str_rcpp(std::string dna_seq, int minSTR = 1, int maxSTR = 6, int minSTRlen = 10, int minReps = 3) {
    int total_bases = dna_seq.length();
    std::vector<REP> srep;
    for (int i = 0; i < (total_bases - minSTRlen); i++) {
        while (i < total_bases - 1 && std::tolower(dna_seq[i]) == 'n') i++;
        for (int rpsz = minSTR; rpsz <= maxSTR; rpsz++) {
            int reps = 1;
            int j = i + rpsz;
            while (j + rpsz <= total_bases && dna_seq.substr(i, rpsz) == dna_seq.substr(j, rpsz)) {
                reps++; j += rpsz;
            }
            if (reps >= minReps) {
                int remainder = 0;
                int rs = i, re = j;
                while (re < total_bases && std::tolower(dna_seq[rs]) == std::tolower(dna_seq[re])) {
                    remainder++; rs++; re++;
                }
                if (((reps * rpsz) + remainder) >= minSTRlen) {
                    if (srep.empty() || srep.back().end < re) {
                        REP r; r.start = i + 1; r.end = re; r.num = reps; r.loop = nonBstr(dna_seq, i + 1, rpsz);
                        r.len = rpsz; r.sub = remainder; r.strand = 0; r.special = 0;
                        srep.push_back(r);
                        i = re - minSTRlen; // Skip ahead
                        break;
                    }
                }
            }
        }
    }
    return reps_to_df(srep);
}

/* Helper for Z-DNA classification */
int pupy(const std::string& dna, int pos) {
    char c1 = std::tolower(dna[pos]);
    char c2 = std::tolower(dna[pos+1]);
    if (c1 == 'a' && c2 == 'c') return 3;
    if (c1 == 't' && c2 == 'g') return 3;
    if (c1 == 'c') {
        if (c2 == 'g') return 25;
        if (c2 == 'a') return 3;
    }
    if (c1 == 'g') {
        if (c2 == 'c') return 25;
        if (c2 == 't') return 3;
    }
    return 0;
}

/* Z-DNA Finder */
// [[Rcpp::export]]
DataFrame find_zdna_rcpp(std::string dna_seq, int minZ = 10) {
    int total_bases = dna_seq.length();
    std::vector<REP> zrep;
    int npy = 1, kvsum = 0;
    int i = 0;
    while (i < (total_bases - 1)) {
        int tmpPPY = pupy(dna_seq, i);
        if (tmpPPY > 0) {
            npy++; kvsum += tmpPPY;
        } else {
            if (npy >= minZ) {
                REP r; r.start = i - npy + 2; r.len = npy; r.loop = kvsum / 2; r.num = 0; r.end = i + 1; r.sub = 0; r.strand = 0; r.special = 0;
                zrep.push_back(r);
            }
            npy = 1; kvsum = 0;
        }
        i++;
    }
    return reps_to_df(zrep);
}

/* Post-processing: remove subsets and overlaps */
void process_repeats_internal(std::vector<REP>& reps, bool centeredOnly = false) {
    if (reps.size() < 2) return;
    std::sort(reps.begin(), reps.end(), [](const REP& a, const REP& b) {
        if (a.start != b.start) return a.start < b.start;
        return b.end < a.end;
    });

    for (size_t i = 1; i < reps.size(); i++) {
        if (reps[i].end <= reps[i-1].end) {
            bool centered = (reps[i].start - reps[i-1].start) == (reps[i-1].end - reps[i].end);
            if (centered) {
                if (reps[i].start == reps[i-1].start) {
                    if (reps[i].loop >= reps[i-1].loop) {
                        reps.erase(reps.begin() + i); i--;
                    } else {
                        reps.erase(reps.begin() + i - 1); i--;
                    }
                } else {
                    reps.erase(reps.begin() + i); i--;
                }
            } else if (!centeredOnly) {
                if (reps[i-1].loop <= reps[i].loop) {
                    reps.erase(reps.begin() + i); i--;
                }
            }
        }
    }
}

/* Classification helper: is_subset logic */
void classify_motifs_internal(std::vector<REP>& reps, const std::string& dna, char type, int max_loop, int limit) {
    for (auto& r : reps) {
        r.special = 0;
        if (type == 'M') { // Mirror repeats -> Triplex
            int nY = 0, nR = 0;
            for (int j = r.start - 1; j < r.start + r.len - 1; j++) {
                char c = std::tolower(dna[j]);
                if (c == 'a' || c == 'g') nR++;
                else if (c == 't' || c == 'c') nY++;
            }
            int Ypercent = (nR == 0) ? 100 : (nY * 100 / (nY + nR));
            if (r.loop < max_loop && Ypercent >= limit) r.special = 1;
        } else if (type == 'D') { // Direct repeats -> Slipped
            if (r.loop <= max_loop) r.special = 1;
        } else if (type == 'Z') { // Z-DNA
            if (r.loop >= limit) r.special = 1;
        } else if (type == 'I') { // Inverted repeats -> Cruciform
            if (r.len >= limit && r.loop <= max_loop) r.special = 1;
        }
    }
}
