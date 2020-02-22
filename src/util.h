#ifndef UTIL_H
#define UTIL_H

void initUses(nle_input_use_t *uses);
void addUses(nle_input_use_t *dest, nle_input_use_t *src);
void printUses(nle_input_use_t *uses);
char *underscore(char *str, int len);
unsigned int gcd(unsigned int u, unsigned int v);
void checkSymmetry2(int *symmetry, int term1, int term2);
void checkSymmetry3(int *symmetry, int term1, int term2, int term3);
int interesting(int range, int max_int, int filter_int, double inld);

#endif // UTIL_H
