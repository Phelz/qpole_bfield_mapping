#include <cmath>

#include "constants.h"


inline void subtract(const PRECISION_TYPE* a, const PRECISION_TYPE* b, PRECISION_TYPE* result) {
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
}

inline void cross_product(const PRECISION_TYPE* a, const PRECISION_TYPE* b, PRECISION_TYPE* result) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

inline PRECISION_TYPE dot_product(const PRECISION_TYPE* a, const PRECISION_TYPE* b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline PRECISION_TYPE norm(const PRECISION_TYPE* a) {
    return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

inline void scale(const PRECISION_TYPE* a, PRECISION_TYPE scalar, PRECISION_TYPE* result) {
    result[0] = scalar * a[0];
    result[1] = scalar * a[1];
    result[2] = scalar * a[2];
}

inline void add(PRECISION_TYPE* a, const PRECISION_TYPE* b) {
    a[0] += b[0];
    a[1] += b[1];
    a[2] += b[2];
}