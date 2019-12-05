#ifndef UTIL_H
#define UTIL_H

void initUses(input_use *uses);
void addUses(input_use *dest, input_use *src);
void printUses(input_use *uses);
char *underscore(char *str, int len);
unsigned int gcd(unsigned int u, unsigned int v);
unsigned int gcd2(unsigned int u, unsigned int v);
void checkSymmetry(int *symmetry, int left, int middle, int right);
int interesting(int range, double inld);

#endif // UTIL_H
