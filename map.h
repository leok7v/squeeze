#ifndef map_header_included
#define map_header_included

#include <stdint.h>

// Braindead hash map for lz77 compression dictionary
// Only supports d from 2 to 255 b because
// 255 to ~2 b compression is about 1% of source and is good enough.

typedef uint8_t map_entry_type[256]; // d[0] number of b [2..127]

typedef struct {
    map_entry_type* entry;
    int32_t n; // entry[n]
    int32_t entries;
    int32_t max_chain;
    int32_t max_bytes;
} map_type;

typedef struct {
    void        (*init)(map_type* m, map_entry_type entry[], size_t n);
    const void* (*data)(const map_type* m, int32_t i);
    uint8_t     (*bytes)(const map_type* m, int32_t i);
    int32_t     (*get)(const map_type* m, const void* data, uint8_t bytes);
    int32_t     (*put)(map_type* m, const void* data, uint8_t bytes);
    int32_t     (*best)(const map_type* m, const void* data, size_t bytes);
    void        (*clear)(map_type *m);
} map_interface;

// map.put()  is no operation if map is filled to 75% or more
// map.get()  returns index of matching entry or -1
// map.best() returns index of longest matching entry or -1

extern map_interface map;

#endif // map_header_included

#if defined(map_implementation) && !defined(map_implemented)

#define map_implemented

#include <string.h>

#ifndef assert
#include <assert.h>
#endif

#ifndef null
#define null ((void*)0)
#endif

// FNV Fowler–Noll–Vo hash function
// https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function

// FNV offset basis for 64-bit:
static const uint64_t map_hash_init = 0xCBF29CE484222325;

// FNV prime for 64-bit
static const uint64_t map_prime64   = 0x100000001B3;

static inline uint64_t map_hash64_byte(uint64_t hash, const uint8_t byte) {
    return (hash ^ (uint64_t)byte) * map_prime64;
}

static inline uint64_t map_hash64(const uint8_t* data, size_t bytes) {
    enum { max_bytes = sizeof(((map_type*)(null))->entry[0]) - 1 };
    assert(2 <= bytes && bytes <= max_bytes);
    uint64_t hash = map_hash_init;
    for (size_t i = 0; i < bytes; i++) {
        hash = map_hash64_byte(hash, data[i]);
    }
    return hash;
}

static void map_init(map_type* m, map_entry_type entry[], size_t n) {
    assert(16 < n && n < 1024 * 1024);
    m->n = (int32_t)n;
    m->entry = entry;
    for (int32_t i = 0; i < m->n; i++) { m->entry[i][0] = 0; }
    m->entries = 0;
    m->max_chain = 0;
    m->max_bytes = 0;
}

static inline const void* map_data(const map_type* m, int32_t i) {
    assert(0 <= i && i < m->n);
    return m->entry[i][0] > 0 ? &m->entry[i][1] : null;
}

static inline uint8_t map_bytes(const map_type* m, int32_t i) {
    assert(0 <= i && i < m->n);
    return m->entry[i][0];
}


static int32_t map_get_hashed(const map_type* m, uint64_t hash,
                              const void* d, uint8_t b) {
    enum { max_bytes = sizeof(m->entry[0]) - 1 };
    assert(2 <= b && b <= max_bytes);
    const map_entry_type* entries = m->entry;
    size_t i = (size_t)hash % m->n;
    // Because map is filled to 3/4 only there will always be
    // an empty slot at the end of the chain.
    while (entries[i][0] > 0) {
        if (entries[i][0] == b && memcmp(entries[i] + 1, d, b) == 0) {
            return (int32_t)i;
        }
        i = (i + 1) % m->n;
    }
    return -1;
}

static int32_t map_get(const map_type* m, const void* d, uint8_t b) {
    return map_get_hashed(m, map_hash64(d, b), d, b);
}

static int32_t map_put(map_type* m, const uint8_t* d, uint8_t b) {
    enum { max_bytes = sizeof(m->entry[0]) - 1 };
    assert(2 <= b && b <= max_bytes);
    if (m->entries < m->n * 3 / 4) {
        map_entry_type* entries = m->entry;
        uint64_t hash = map_hash64(d, b);
        size_t i = (size_t)hash % m->n;
        int32_t chain = 0; // max chain length
        while (entries[i][0] > 0) {
            if (entries[i][0] == b && memcmp(entries[i] + 1, d, b) == 0) {
                return (int32_t)i; // found match with existing entry
            }
            chain++;
            i = (i + 1) % m->n;
            assert(chain < m->n); // looping endlessly?
        }
        if (chain > m->max_chain) { m->max_chain = chain; }
        if (b  > m->max_bytes) { m->max_bytes = b; }
        entries[i][0] = b;
        memcpy(entries[i] + 1, d, b);
        m->entries++;
        return (int32_t)i;
    }
    return -1;
}

static int32_t map_best(const map_type* m, const void* data, size_t bytes) {
    enum { max_bytes = sizeof(m->entry[0]) - 1 };
    int32_t best = -1; // best (longest) result
    if (bytes > 1) {
        const uint8_t  b = (uint8_t)(bytes <= max_bytes ? bytes : max_bytes);
        const uint8_t* d = (uint8_t*)data;
        uint64_t hash = map_hash64_byte(map_hash_init, d[0]);
        for (uint8_t i = 1; i < b - 1; i++) {
            hash = map_hash64_byte(hash, d[i]);
            int32_t r = map_get_hashed(m, hash, data, i + 1);
            if (r != -1) {
                best = r;
            } else if (best != -1) {
                break; // will return longest matching entry index
            }
        }
    }
    return best;
}

static void map_clear(map_type *m) {
    for (int32_t i = 0; i < m->n; i++) {
        m->entry[i][0] = 0;
    }
    m->entries = 0;
    m->max_chain = 0;
    m->max_bytes = 0;
}

map_interface map = {
    .init  = map_init,
    .data  = map_data,
    .bytes = map_bytes,
    .get   = map_get,
    .put   = map_put,
    .best  = map_best,
    .clear = map_clear
};

#endif // map_implementation

