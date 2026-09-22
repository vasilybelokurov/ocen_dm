// POSIX aligned allocation for falcON headers on ARM GCC, which lacks this
// x86 convenience header. No SIMD operations or changes to falcON's ABI.
#ifndef OCEN_MM_MALLOC_H
#define OCEN_MM_MALLOC_H
#include <cstdlib>
static inline void* _mm_malloc(std::size_t size, std::size_t alignment) {
    if (alignment <= sizeof(void*)) alignment = sizeof(void*);
    void* result = 0;
    return posix_memalign(&result, alignment, size) == 0 ? result : 0;
}
static inline void _mm_free(void* pointer) { std::free(pointer); }
#endif
