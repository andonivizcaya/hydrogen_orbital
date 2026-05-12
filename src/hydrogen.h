// a lot of stuff was stolen from musl
#ifndef HYDROGEN_H
#define HYDROGEN_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    double *items;
    size_t count;
    size_t capacity;
} HydrogenRow;

typedef struct {
    HydrogenRow *items;
    size_t count;
    size_t capacity;
} HydrogenMatrix;

#define DA_INIT_CAP 256

#ifdef HYDROLIB_WEB

extern unsigned char __heap_base;

#ifndef WASM_PAGE_SIZE
#define WASM_PAGE_SIZE 65536u
#endif // WASM_PAGE_SIZE

static inline size_t wasm_mem_byte_size(void)
{
    return (size_t)__builtin_wasm_memory_size(0u)*(size_t)WASM_PAGE_SIZE;
}

static inline int wasm_mem_grow_to(size_t addr_end)
{
    size_t mem_end = wasm_mem_byte_size();

    if (addr_end <= mem_end) return 0;

    size_t deficit          = addr_end - mem_end;
    uint32_t pages          = (uint32_t)((deficit + WASM_PAGE_SIZE - 1u)/WASM_PAGE_SIZE);
    unsigned long long prev = __builtin_wasm_memory_grow(0u, pages);

    return prev != (unsigned long long)-1 ? 0 : -1;
}

static size_t g_bump;

static inline void *header_ptr(void *p)
{
    return (char *)p - sizeof(size_t);
}

void *wasm_malloc(size_t n)
{
    if (n == 0) n = 1;

    const size_t align = sizeof(size_t);
    uintptr_t hb       = (uintptr_t)&__heap_base;
    uintptr_t raw      = hb + g_bump;

    uintptr_t user       = (raw + sizeof(size_t) + align - 1u) & ~(uintptr_t)(align - 1u);
    uintptr_t block_end  = user + n;

    if (wasm_mem_grow_to(block_end) != 0) return NULL;

    *(size_t *)(user - sizeof(size_t)) = n;
    g_bump                             = (size_t)(block_end - hb);

    return (void *)user;
}

void wasm_free(void *p)
{
    (void)p;
}

static inline void *memcpy(void *restrict dest, const void *restrict src, size_t n)
{
    unsigned char *d       = dest;
    const unsigned char *s = src;
    for (; n; n--) *d++ = *s++;
    return dest;
}

void *wasm_realloc(void *old, size_t new_nbytes)
{
    if (new_nbytes == 0) {
        wasm_free(old);
        return NULL;
    }
    if (!old) return wasm_malloc(new_nbytes);

    size_t old_nbytes = *((size_t *)header_ptr(old));

    if (new_nbytes <= old_nbytes) {
        *((size_t *)header_ptr(old)) = new_nbytes;
        return old;
    }

    void *fresh = wasm_malloc(new_nbytes);
    if (!fresh) return NULL;

    unsigned char *d = fresh;
    unsigned char *s = old;

    for (size_t i = 0; i < old_nbytes; i++) {
        d[i] = s[i];
    }

    wasm_free(old);
    return fresh;
}

static inline size_t hydrogen_cstrlen(const char *s)
{
    size_t n = 0;
    while (s[n]) n++;
    return n;
}

#define ASSERT(cond)          \
    do {                      \
        if (!(cond)) {        \
            __builtin_trap(); \
        }                     \
    } while (0)

#define FREE      wasm_free
#define MALLOC    wasm_malloc
#define MEMCPY(d, s, n) memcpy((d), (s), (n))
#define REALLOC   wasm_realloc
#define STRLEN(s) hydrogen_cstrlen((s))

__attribute__((noreturn)) static inline void exit(int status)
{
    (void)status;
    __builtin_trap();
}

#define sinf(x) __builtin_sinf((float)(x))
#define cosf(x) __builtin_cosf((float)(x))
#define tanf(x) __builtin_tanf((float)(x))
#define sqrtf(x) __builtin_sqrtf((float)(x))
#define fabsf(x) __builtin_fabsf((float)(x))
#define floorf(x) __builtin_floorf((float)(x))
#define fminf(a, b) __builtin_fminf((float)(a), (float)(b))
#define fmaxf(a, b) __builtin_fmaxf((float)(a), (float)(b))
#define atan2f(y, x) __builtin_atan2f((float)(y), (float)(x))
#define acosf(x) __builtin_acosf((float)(x))
#define asinf(x) __builtin_asinf((float)(x))
#define roundf(x) __builtin_roundf((float)(x))
#define sin(x) __builtin_sin((double)(x))
#define cos(x) __builtin_cos((double)(x))
#define tan(x) __builtin_tan((double)(x))
#define sqrt(x) __builtin_sqrt((double)(x))
#define fabs(x) __builtin_fabs((double)(x))
#define pow(x, y) __builtin_pow((double)(x), (double)(y))
#define powf(x, y) __builtin_powf((float)(x), (float)(y))
#define fmax(a, b) __builtin_fmax((double)(a), (double)(b))
#define fmin(a, b) __builtin_fmin((double)(a), (double)(b))

#else
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#define ASSERT assert
#define FREE free
#define MALLOC malloc
#define MEMCPY(d, s, n) memcpy((d), (s), (n))
#define REALLOC realloc
#define STRLEN(s) strlen((s))
#endif // HYDROLIB_WEB

static inline int hydrogen_iabs(int x)
{
    return x < 0 ? -x : x;
}

#define max(a, b) a > b ? a : b

#define da_reserve(da, expected_capacity)                                              \
    do {                                                                               \
        if ((expected_capacity) > (da)->capacity) {                                    \
            if ((da)->capacity == 0) {                                                 \
                (da)->capacity = DA_INIT_CAP;                                          \
            }                                                                          \
            while ((expected_capacity) > (da)->capacity) {                             \
                (da)->capacity *= 2;                                                  \
            }                                                                          \
            (da)->items = REALLOC((da)->items, (da)->capacity * sizeof(*(da)->items)); \
            ASSERT((da)->items != NULL && "Buy more RAM lol");                         \
        }                                                                              \
    } while (0)

#define da_append(da, item)                  \
    do {                                     \
        da_reserve((da), (da)->count + 1);   \
        (da)->items[(da)->count++] = (item); \
    } while (0)

#define da_free(da) FREE((da).items)
#define da_foreach(Type, it, da) for (Type *it = (da)->items; it < (da)->items + (da)->count; ++it)
#endif // HYDROGEN_H


#if defined(HYDROLIB_WEB) && defined(RAYLIB_H) && defined(RAYMATH_H)

#ifndef HYDROGEN_RWEB
#define HYDROGEN_RWEB

void *MemAlloc(unsigned int size)
{
    return MALLOC((size_t)size);
}

void *MemRealloc(void *ptr, unsigned int size)
{
    return REALLOC(ptr, (size_t)size);
}

void MemFree(void *ptr)
{
    FREE(ptr);
}

#ifndef RL_CULL_DISTANCE_NEAR
#define RL_CULL_DISTANCE_NEAR 0.05f
#endif // RL_CULL_DISTANCE_NEAR

#ifndef RL_CULL_DISTANCE_FAR
#define RL_CULL_DISTANCE_FAR 4000.0f
#endif // RL_CULL_DISTANCE_FAR

#define RCAMERA_IMPLEMENTATION
#include "rcamera.h"

__attribute__((import_module("env"), import_name("hydrolib_js_set_entry")))
extern void hydrolib_js_set_entry(void (*fn)(void));

#endif // HYDROGEN_RWEB

#endif // defined(HYDROLIB_WEB) && defined(RAYLIB_H) && defined(RAYMATH_H)