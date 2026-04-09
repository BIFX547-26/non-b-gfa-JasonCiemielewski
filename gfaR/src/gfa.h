#ifndef GFA_H_
#define GFA_H_

#include <vector>
#include <string>

/* Structure definition for REP */
typedef struct REP {
    int start;
    int loop;
    int len;
    int num;
    int end;
    int sub;
    int strand; // 0 for plus, 1 for minus
    int special; // 1 for yes, 0 for no
} REP;

/* A-Tract */
typedef struct A_Tract {
    int strt;
    short len;
} A_Tract;

/* potential Bent DNA */
typedef struct potential_Bent_DNA {
    double a_center;
    int strt;
    int end;
} potential_Bent_DNA;

/* G-Island */
typedef struct G_Island {
    int strt;
    int len;
} G_Island;

/* Utility functions */
std::string complement(const std::string& dna);
std::string reverse_complement(const std::string& dna);

#endif /* GFA_H_ */
