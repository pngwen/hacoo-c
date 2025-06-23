#ifndef VECTOR_H
#define VECTOR_H

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define VECTOR_INIT_CAPACITY 4

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
    vec.data = malloc(VECTOR_INIT_CAPACITY * sizeof(TYPE));                    \
    vec.size = 0;                                                               \
    vec.capacity = VECTOR_INIT_CAPACITY;                                        \
    return vec;                                                                 \
}                                                                               \
                                                                                \
static inline void NAME##_push_back(NAME *vec, TYPE value) {                    \
    if (vec->size == vec->capacity) {                                           \
        vec->capacity *= 2;                                                     \
        vec->data = realloc(vec->data, vec->capacity * sizeof(TYPE));          \
        assert(vec->data != NULL);                                              \
    }                                                                           \
    vec->data[vec->size++] = value;                                             \
}                                                                               \
                                                                                \
static inline TYPE NAME##_get(NAME *vec, size_t index) {                        \
    assert(index < vec->size);                                                  \
    return vec->data[index];                                                    \
}                                                                               \
                                                                                \
static inline void NAME##_free(NAME *vec) {                                     \
    free(vec->data);                                                            \
    vec->data = NULL;                                                           \
    vec->size = vec->capacity = 0;                                              \
}
#endif // VECTOR_H
