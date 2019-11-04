#include "cdstd.h"
#include <immintrin.h>
using namespace std;
__m256d sub(__m256d x)
{
    return _mm256_mul_pd(x, x);
}
