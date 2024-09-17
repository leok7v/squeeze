#ifndef squeeze_header_included
#define squeeze_header_included

#include "bitstream.h"
#include "deflate.h"
#include "huffman.h"
#include "map.h"

enum { // "nyt" stands for Not Yet Transmitted (see Vitter Algorithm)
    squeeze_min_win_bits  =  10,
    squeeze_max_win_bits  =  15,
    squeeze_min_map_bits  =  16,
    squeeze_max_map_bits  =  28,
    squeeze_lit_nyt       =  deflate_sym_max + 1,
    squeeze_pos_nyt       =  deflate_pos_max + 1,
};

typedef struct {
    errno_t error; // sticky
    huffman_tree_type  lit; // 0..255 literal bytes; 257-285 length
    huffman_tree_type  pos; // positions tree of 1^win_bits
    huffman_node_type* lit_nodes; // 512 * 2 - 1
    huffman_node_type* pos_nodes; //  32 * 2 - 1 
    map_type           map;
    map_entry_type*    map_entry;
    bitstream_type*    bs;
} squeeze_type;

#define squeeze_size_mul(name, count) (                                \
    ((uint64_t)(count) >= ((SIZE_MAX / 4) / (uint64_t)sizeof(name))) ? \
    0 : (size_t)((uint64_t)sizeof(name) * (uint64_t)(count))           \
)

#define squeeze_size_implementation(map_bits) (                        \
    (sizeof(squeeze_type)) +                                           \
    /* lit_nodes: */                                                   \
    squeeze_size_mul(huffman_node_type, (512uLL * 2ULL - 1ULL)) +      \
    /* pos_nodes: */                                                   \
    squeeze_size_mul(huffman_node_type, ((1uLL << 5) * 2ULL - 1ULL)) + \
    /* pos_nodes: */                                                   \
    squeeze_size_mul(map_entry_type, (1uLL << (map_bits)))             \
)

#define squeeze_sizeof(win_bits, map_bits) (                           \
    (sizeof(size_t) == sizeof(uint64_t)) &&                            \
    (squeeze_min_win_bits <= (win_bits)) &&                            \
                             ((win_bits) <= squeeze_max_win_bits) ?    \
    (size_t)squeeze_size_implementation((map_bits)) : 0                \
)

typedef struct {
    squeeze_type* (*alloc)(bitstream_type* bs, uint8_t win_bits,
                           uint8_t map_bits);
    errno_t (*init_with)(squeeze_type* s, void* memory, size_t size,
                         uint8_t win_bits, uint8_t map_bits);
    // `win_bits` is a log2 of window size in bytes in range
    // [squeeze_min_win_bits..squeeze_max_win_bits]
    void (*write_header)(bitstream_type* bs, uint64_t bytes, uint8_t win_bits);
    void (*compress)(squeeze_type* s, const uint8_t* data, size_t bytes, 
                     uint16_t window);
    void (*read_header)(bitstream_type* bs, uint64_t *bytes, uint8_t *win_bits);
    void (*decompress)(squeeze_type* s, uint8_t* data, size_t bytes);
    void (*free)(squeeze_type* s);
} squeeze_interface;

extern squeeze_interface squeeze;

#endif // squeeze_header_included

#if defined(squeeze_implementation) && !defined(squeeze_implemented)

#define squeeze_implemented

#define bitstream_implementation
#include "bitstream.h"

#define huffman_implementation
#include "huffman.h"

#ifndef null
#define null ((void*)0) // like null_ptr better than NULL (0)
#endif

#ifndef countof
#define countof(a) (sizeof(a) / sizeof((a)[0]))
#endif

#ifndef assert
#include <assert.h>
#endif

static squeeze_type* squeeze_alloc(bitstream_type* bs, uint8_t win_bits,
                                   uint8_t map_bits) {
    const uint64_t bytes = squeeze_sizeof(win_bits, map_bits);
    squeeze_type* s = (squeeze_type*)calloc(1, (size_t)bytes);
    if (s != null) {
        squeeze.init_with(s, s, bytes, win_bits, map_bits);
        s->bs = bs;
    }
    return s;
}

static void squeeze_free(squeeze_type* s) {
    free(s);
}

static errno_t squeeze_init_with(squeeze_type* s, void* memory, size_t size,
                                 uint8_t win_bits, uint8_t map_bits) {
    errno_t r = 0;
    assert(squeeze_min_win_bits <= win_bits && win_bits <= squeeze_max_win_bits);
    assert(map_bits == 0 ||
           squeeze_min_map_bits <= map_bits && map_bits <= squeeze_max_map_bits);
    size_t expected = squeeze_sizeof(win_bits, map_bits);
    // 167,936,192 bytes for (win_bits = 11, map_bits = 19)
    assert(size == expected);
    if (expected == 0 || memory == null || size != expected) {
        r = EINVAL;
    } else {
        uint8_t* p = (uint8_t*)memory;
        memset(memory, 0, sizeof(squeeze_type));
        p += sizeof(squeeze_type);
        const size_t lit_n = 512; // literals in range [0..286] always 512
        const size_t pos_n = 32;
        const size_t map_n = 1uLL << map_bits;
        const size_t lit_m = lit_n * 2 - 1;
        const size_t pos_m = pos_n * 2 - 1;
        s->lit_nodes = (huffman_node_type*)p; p += sizeof(huffman_node_type) * lit_m;
        s->pos_nodes = (huffman_node_type*)p; p += sizeof(huffman_node_type) * pos_m;
        s->map_entry = (map_entry_type*)p;    p += sizeof(map_entry_type)    * map_n;
        assert(p == (uint8_t*)memory + size);
        huffman_init(&s->lit, s->lit_nodes, lit_m);
        huffman_init(&s->pos, s->pos_nodes, pos_m);
        if (map_bits != 0) {
            map_init(&s->map, s->map_entry, map_n);
        } else { // map is not used in decompress
            memset(&s->map, 0x00, sizeof(s->map));
        }
    }
    return r;
}

#define squeeze_if_error_return(s) do { \
    if (s->error) { return; }           \
} while (0)

#define squeeze_return_invalid(s) do {  \
    s->error = EINVAL;                  \
    return;                             \
} while (0)

static inline void squeeze_write_bit(squeeze_type* s, bool bit) {
    if (s->error == 0) {
        bitstream_write_bit(s->bs, bit);
        s->error = s->bs->error;
    }
}

static inline void squeeze_write_bits(squeeze_type* s,
                                      uint64_t b64, uint8_t bits) {
    if (s->error == 0) {
        bitstream_write_bits(s->bs, b64, bits);
        s->error = s->bs->error;
    }
}

static inline void squeeze_write_huffman(squeeze_type* s, huffman_tree_type* t,
                                         int32_t i) {
    assert(t != null && t->node != null);
    assert(0 <= i && i < t->n); // leaf symbol (literal)
    assert(1 <= t->node[i].bits && t->node[i].bits < 64);
    squeeze_write_bits(s, t->node[i].path, (uint8_t)t->node[i].bits);
    huffman_inc_frequency(t, i); // after the path is written
}

static inline void squeeze_flush(squeeze_type* s) {
    if (s->error == 0) {
        bitstream_flush(s->bs);
        s->error = s->bs->error;
    }
}

static void squeeze_write_header(bitstream_type* bs, uint64_t bytes,
                                 uint8_t win_bits) {
    if (win_bits < squeeze_min_win_bits || win_bits > squeeze_max_win_bits) {
        bs->error = EINVAL;
    } else {
        enum { bits64 = sizeof(uint64_t) * 8 };
        bitstream_write_bits(bs, (uint64_t)bytes, bits64);
        enum { bits8 = sizeof(uint8_t) * 8 };
        bitstream_write_bits(bs, win_bits, bits8);
    }
}

static uint8_t squeeze_log2_of_pow2(uint64_t pow2) {
    assert(pow2 > 0 && (pow2 & (pow2 - 1)) == 0);
    if (pow2 > 0 && (pow2 & (pow2 - 1)) == 0) {
        uint8_t bit = 0;
        while (pow2 >>= 1) { bit++; }
        return bit;
    } else {
        return 0xFF; // error
    }
}

static inline void squeeze_encode_literal(squeeze_type* s, const uint16_t lit) {
    assert(0 <= lit && lit <= 285);
    if (s->lit.node[lit].bits == 0) {
        assert(s->lit.node[squeeze_lit_nyt].bits != 0);
        squeeze_write_huffman(s, &s->lit, squeeze_lit_nyt);
        squeeze_write_bits(s, lit, 9);
        if (!huffman_insert(&s->lit, lit)) { s->error = E2BIG; }
    } else {
        squeeze_write_huffman(s, &s->lit, lit);
    }
}

static inline void squeeze_encode_len(squeeze_type* s, const uint16_t len) {
    assert(3 <= len && len < deflate_len_count);
    uint8_t  i = deflate_len_index[len];
    uint16_t b = deflate_len_base[i];
    uint8_t  x = deflate_len_extra_bits[i];
    squeeze_encode_literal(s, (uint16_t)(257 + i));
    assert(b <= len && len - b <= (uint16_t)(1u << x));
    if (x > 0) { squeeze_write_bits(s, (len - b), x); }
}

static inline void squeeze_encode_pos(squeeze_type* s, const uint16_t pos) {
    assert(0 <= pos && pos <= 0x7FFF);
    uint8_t  i = deflate_pos_index[pos];
    uint16_t b = deflate_pos_base[i];
    uint8_t  x = deflate_pos_extra_bits[i];
    if (s->pos.node[i].bits == 0) {
        assert(s->pos.node[squeeze_pos_nyt].bits != 0);
        squeeze_write_huffman(s, &s->pos, squeeze_pos_nyt);
        squeeze_write_bits(s, i, 5); // 0..29
        if (!huffman_insert(&s->pos, i)) { s->error = E2BIG; }
    } else {
        squeeze_write_huffman(s, &s->pos, i);
    }
    assert(b <= pos && pos - b <= (uint16_t)(1u << x));
    if (x > 0) { squeeze_write_bits(s, (pos - b), x); }
}

// #define SQUEEZE_MAP_STATS // uncomment to print stats

static void squeeze_compress(squeeze_type* s, const uint8_t* data, uint64_t bytes, 
                             uint16_t window) {
    #ifdef SQUEEZE_MAP_STATS
        static double   map_distance_sum;
        static double   map_len_sum;
        static uint64_t map_count;
        map_distance_sum = 0;
        map_len_sum = 0;
        map_count = 0;
        size_t br_bytes = 0; // source bytes encoded as back references
        size_t li_bytes = 0; // source bytes encoded "as is" literals
    #endif
    huffman_insert(&s->lit, squeeze_lit_nyt);
    huffman_insert(&s->pos, squeeze_pos_nyt);
    deflate_init();
    size_t i = 0;
    while (i < bytes && s->error == 0) {
        size_t len = 0;
        size_t pos = 0;
        if (i >= 1) {
            size_t j = i - 1;
            size_t min_j = i >= window ? i - window + 1 : 0;
            for (;;) {
                assert((i - j) < window);
                const size_t n = bytes - i;
                size_t k = 0;
                while (k < n && data[j + k] == data[i + k] && k < deflate_len_max) {
                    k++;
                }
                if (k >= deflate_len_min && k > len) {
                    len = k;
                    pos = i - j;
                    if (len == deflate_len_max) { break; }
                }
                if (j == min_j) { break; }
                j--;
            }
        }
        if (s->map.n > 0) {
            int32_t  best = map_best(&s->map, data + i, bytes - i);
            uint32_t best_bytes = 
                     best < 0 ? 0 : s->map.entry[best].bytes;
            uint32_t best_distance = best < 0 ? 
                     0 : (uint32_t)(data + i - s->map.entry[best].data);
            if (best_distance < 0x7FFF && best_bytes > len && best_bytes > 4) {
                assert(best_bytes >= deflate_len_min);
                assert(memcmp(data + i - best_distance, data + i, best_bytes) == 0);
                len = best_bytes;
                pos = best_distance;
                #ifdef SQUEEZE_MAP_STATS
                    map_distance_sum += pos;
                    map_len_sum += len;
                    map_count++;
                #endif
            }
        }
        if (len >= deflate_len_min) {
            assert(0 < pos && pos <= 0x7FFF);
            squeeze_encode_len(s, (uint16_t)len);
            squeeze_encode_pos(s, (uint16_t)pos);
            if (s->map.n > 0) {
                map_put(&s->map, data + i, (uint32_t)len);
            }
            i += len;
            #ifdef SQUEEZE_MAP_STATS
                br_bytes += len;
            #endif
        } else {
            #ifdef SQUEEZE_MAP_STATS
            li_bytes++;
            #endif
            squeeze_encode_literal(s, data[i]);
            i++;
        }
    }
    squeeze_flush(s);
    #ifdef SQUEEZE_MAP_STATS
        double br_percent = (100.0 * br_bytes) / (br_bytes + li_bytes);
        double li_percent = (100.0 * li_bytes) / (br_bytes + li_bytes);
        printf("entropy literals: %.2f %.2f%% back references: %.2f %.2f%%\n",
            huffman_entropy(&s->lit), li_percent,
            huffman_entropy(&s->pos), br_percent);
        if (map_count > 0) {
            printf("avg dic distance: %.1f length: %.1f count: %lld\n", 
                    map_distance_sum / map_count,
                    map_len_sum / map_count, map_count);
        }
    #endif
}

static inline uint64_t squeeze_read_bit(squeeze_type* s) {
    bool bit = 0;
    if (s->error == 0) {
        bit = bitstream_read_bit(s->bs);
        s->error = s->bs->error;
    }
    return bit;
}

static inline uint64_t squeeze_read_bits(squeeze_type* s, uint32_t n) {
    assert(n <= 64);
    uint64_t bits = 0;
    if (s->error == 0) {
        bits = bitstream_read_bits(s->bs, n);
    }
    return bits;
}

static inline uint64_t squeeze_read_huffman(squeeze_type* s, huffman_tree_type* t) {
    const int32_t m = t->n * 2 - 1;
    int32_t i = m - 1; // root
    bool bit = squeeze_read_bit(s);
    while (s->error == 0) {
        i = bit ? t->node[i].rix : t->node[i].lix;
        assert(0 <= i && i < m);
        if (t->node[i].lix < 0 && t->node[i].rix < 0) { break; } // leaf
        bit = squeeze_read_bit(s);
    }
    assert(0 <= i && i < t->n); // leaf symbol (literal)
    huffman_inc_frequency(t, i);
    return (uint64_t)i;
}

static void squeeze_read_header(bitstream_type* bs, uint64_t *bytes,
                                uint8_t *win_bits) {
    uint64_t b = bitstream_read_bits(bs, sizeof(uint64_t) * 8);
    uint64_t w = bitstream_read_bits(bs, sizeof(uint8_t) * 8);
    if (bs->error == 0) {
        if (w < squeeze_min_win_bits || w > squeeze_max_win_bits) {
            bs->error = EINVAL;
        } else if (bs->error == 0) {
            *bytes = b;
            *win_bits = (uint8_t)w;
        }
    }
}

static uint16_t squeeze_read_length(squeeze_type* s, uint16_t lit) {
    const uint8_t base = (uint8_t)(lit - 257);
    if (base >= countof(deflate_len_base)) { 
        s->error = EINVAL;
        return 0;
    } else {
        const uint8_t bits = deflate_len_extra_bits[base];
        if (bits != 0) {
            uint64_t extra = squeeze_read_bits(s, bits);
            assert(deflate_len_base[base] + extra <= UINT16_MAX);
            return s->error == 0 ? 
                (uint16_t)deflate_len_base[base] + (uint16_t)extra : 0;
        } else {
            return deflate_len_base[base];
        }
    }
}

static uint32_t squeeze_read_pos(squeeze_type* s) {
    uint32_t pos = 0;
    uint64_t base = squeeze_read_huffman(s, &s->pos);
    if (s->error == 0 && base == squeeze_pos_nyt) {
        base = squeeze_read_bits(s, 5);
        if (s->error == 0) { 
            if (!huffman_insert(&s->pos, (int32_t)base)) {
                s->error = E2BIG;
            }
        }
    }
    if (s->error == 0 && base >= countof(deflate_pos_base)) {
        s->error = EINVAL;
    } else {
        pos = deflate_pos_base[base];
    }
    if (s->error == 0) { 
        uint8_t bits  = deflate_pos_extra_bits[base];
        if (bits > 0) {
            uint64_t extra = squeeze_read_bits(s, bits);
            if (s->error == 0) { pos += (uint32_t)extra; }
        }
    }
    return pos;
}

static void squeeze_decompress(squeeze_type* s, uint8_t* data, uint64_t bytes) {
    huffman_insert(&s->lit, squeeze_lit_nyt);
    huffman_insert(&s->pos, squeeze_pos_nyt);
    deflate_init();
    size_t i = 0; // output b64[i]
    while (i < bytes && s->error == 0) {
        uint64_t lit = squeeze_read_huffman(s, &s->lit);
        if (s->error != 0) { break; }
        if (lit == squeeze_lit_nyt) {
            lit = squeeze_read_bits(s, 9);
            if (s->error != 0) { break; }
            // TODO: error if insert fails
            if (!huffman_insert(&s->lit, (int32_t)lit)) {
                s->error = E2BIG;
                break;
            }
        }
        if (lit <= 0xFF) {
            data[i] = (uint8_t)lit;
            i++;
        } else {
            assert(257 <= lit && lit < squeeze_lit_nyt);
            if (257 <= lit && lit <= squeeze_lit_nyt) {
                uint32_t len = squeeze_read_length(s, (uint16_t)lit);
                if (s->error != 0) { break; }
                assert(deflate_len_min <= len && len <= deflate_len_max); 
                if (deflate_len_min <= len && len <= deflate_len_max) { 
                    uint32_t pos = squeeze_read_pos(s);
                    if (s->error != 0) { break; }
                    assert(0 < pos && pos <= 0x7FFF);
                    if (0 < pos && pos <= 0x7FFF) {
                        // Cannot do memcpy() because of overlapped regions.
                        // memcpy() may read more than one byte at a time.
                        uint8_t* d = data - (size_t)pos;
                        const size_t n = i + (size_t)len;
                        while (i < n) { data[i] = d[i]; i++; }
                    } else {
                        s->error = EINVAL;
                    }
                } else {
                    s->error = EINVAL; 
                }
            } else {
                s->error = EINVAL;
            }
        }
    }
}

squeeze_interface squeeze = {
    .alloc        = squeeze_alloc,
    .init_with    = squeeze_init_with,
    .write_header = squeeze_write_header,
    .compress     = squeeze_compress,
    .read_header  = squeeze_read_header,
    .decompress   = squeeze_decompress,
    .free         = squeeze_free
};

#endif // squeeze_implementation
