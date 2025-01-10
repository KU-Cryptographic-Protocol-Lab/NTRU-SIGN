#include <stdint.h>
#include "params.h"
#include "reduce.h"
#include <stdlib.h>


int16_t montgomery_reduce(int32_t a)
{
  int16_t t;

  t = (int16_t)a*QINV;
  t = (a - (int32_t)t*NTRUPLUS_Q) >> 16;
  return t;
}

int16_t barrett_reduce(int16_t a)
{
  int16_t t;
  const int16_t v = ((1<<26) + NTRUPLUS_Q/2)/NTRUPLUS_Q;

  t  = ((int32_t)v*a + (1<<25)) >> 26;
  t *= NTRUPLUS_Q;
  return a - t;
}

int16_t caddq(int16_t a) {
    a += (a >> 15) & NTRUPLUS_Q;
    return a;
}

int16_t csubq(int16_t a) {
    int16_t t = ((a-NTRUPLUS_Q) >> 15);
    a -= ~t & NTRUPLUS_Q;
    return a;
}

int16_t caddp(int16_t a) {
    a += (a >> 15) & NTRUPLUS_P;
    return a;
}

int16_t csubp(int16_t a) {
    int16_t t = ((a-NTRUPLUS_P) >> 15);
    a -= ~t & NTRUPLUS_P;
    return a;
}

int max(int16_t a, int16_t b) {
  int16_t t = (a - b) >> 15;
  return (a & ~t) | (b & t);
}

int16_t freezep(int16_t a) {
    a += ((a + halfp) >> 15) & NTRUPLUS_P;
    a -= ((halfp - a) >> 15) & NTRUPLUS_P;
    return a;
}



