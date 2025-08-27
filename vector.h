#ifndef VECTOR_H
#define VECTOR_H

#include "common.cpp"  // brings in MALLOC / FREE / alignment helpers
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define VECTOR_INIT_CAPACITY 4

static inline void* AlignedRealloc(void* ptr, size_t old_size, size_t new_size) {
    void* new_ptr = MALLOC(new_size);
    if (!new_ptr) return NULL;
    if (ptr) {
        memcpy(new_ptr, ptr, old_size < new_size ? old_size : new_size);
        FREE(ptr);
    }
    return new_ptr;
}

// Define a vector for a specific type
#define DEFINE_VECTOR_TYPE(TYPE, NAME)                                          \
typedef struct {                                                                \
    TYPE *data;                                                                 \
    size_t size;                                                                \
    size_t capacity;                                                            \
} NAME;                                                                         \
                                                                                \
static inline NAME NAME##_create() {                                            \
    NAME vec;                                                                   \
    vec.data = (TYPE*) MALLOC(VECTOR_INIT_CAPACITY * sizeof(TYPE));             \
    assert(vec.data != NULL);                                                   \
    vec.size = 0;                                                               \
    vec.capacity = VECTOR_INIT_CAPACITY;                                        \
    return vec;                                                                 \
}                                                                               \
                                                                                \
static inline void NAME##_push_back(NAME *vec, TYPE value) {                    \
    if (vec->size == vec->capacity) {                                           \
        size_t old_size = vec->capacity * sizeof(TYPE);                         \
        vec->capacity *= 2;                                                     \
        vec->data = (TYPE*) AlignedRealloc(vec->data, old_size,                 \
                                   vec->capacity * sizeof(TYPE));               \
        assert(vec->data != NULL);                                              \
    }                                                                           \
    vec->data[vec->size++] = value;                                             \
}                                                                               \
                                                                                \
static inline TYPE NAME##_get(const NAME *vec, size_t index) {                  \
    assert(index < vec->size);                                                  \
    return vec->data[index];                                                    \
}                                                                               \
                                                                                \
static inline void NAME##_free(NAME *vec) {                                     \
    FREE(vec->data);                                                            \
    vec->data = NULL;                                                           \
    vec->size = 0;                                                              \
    vec->capacity = 0;                                                          \
}

#endif // VECTOR_H