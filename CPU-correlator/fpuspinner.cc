
#include <emmintrin.h>

void* calcMaxFlops(void* data)
{
  __m128 a = _mm_set_ps1(1.0);
  __m128 b = _mm_set_ps1(1.0);
  __m128 c = _mm_set_ps1(1.0);
  __m128 d = _mm_set_ps1(1.0);
  __m128 e = _mm_set_ps1(1.0);
  __m128 f = _mm_set_ps1(1.0);
  __m128 g = _mm_set_ps1(1.0);
  __m128 h = _mm_set_ps1(1.0);
  __m128 i = _mm_set_ps1(0.0);
  __m128 j = _mm_set_ps1(0.0);
  __m128 k = _mm_set_ps1(0.0);
  __m128 l = _mm_set_ps1(0.0);
  __m128 m = _mm_set_ps1(0.0);
  __m128 n = _mm_set_ps1(0.0);
  __m128 o = _mm_set_ps1(0.0);
  __m128 p = _mm_set_ps1(0.0);

  volatile __m128 result;

  for (unsigned long long x = 0; x < 1000000000L; x ++) {
    a *= a;
    b *= b;
    c *= c;
    d *= d;
    e *= e;
    f *= f;
    g *= g;
    h *= h;
    i += i;
    j += j;
    k += k;
    l += l;
    m += m;
    n += n;
    o += o;
    p += p;
  }

  result = a + b + c + d + e + f + g + h + i + j + k + l + m + n + o + p;

  return 0;
}
