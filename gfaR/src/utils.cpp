#include "gfa.h"
#include <algorithm>

std::string complement(const std::string& dna) {
    std::string res = dna;
    for (size_t i = 0; i < dna.length(); ++i) {
        switch (dna[i]) {
            case 'a': res[i] = 't'; break;
            case 'c': res[i] = 'g'; break;
            case 'g': res[i] = 'c'; break;
            case 't': res[i] = 'a'; break;
            case 'n': res[i] = 'n'; break;
            case 'A': res[i] = 'T'; break;
            case 'C': res[i] = 'G'; break;
            case 'G': res[i] = 'C'; break;
            case 'T': res[i] = 'A'; break;
            case 'N': res[i] = 'N'; break;
            default: res[i] = dna[i]; break;
        }
    }
    return res;
}

std::string reverse_complement(const std::string& dna) {
    std::string res = complement(dna);
    std::reverse(res.begin(), res.end());
    return res;
}
