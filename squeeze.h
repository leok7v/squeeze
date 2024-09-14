#ifndef squeeze_header_included
#define squeeze_header_included

#include "bitstream.h"
#include "huffman.h"

enum {
    squeeze_min_win_bits  =  10,
    squeeze_max_win_bits  =  15,
    squeeze_lit_nyt       = 286, // Not Yet Transmitted (see Vitter Algorithm)
    squeeze_pos_nyt       =  31  // Not Yet Transmitted
};

typedef struct {
    errno_t error; // sticky
    huffman_tree_type  lit; // 0..255 litteral bytes; 257-285 length
    huffman_tree_type  pos; // positions tree of 1^win_bits
    huffman_node_type* lit_nodes; // 512 * 2 - 1
    huffman_node_type* pos_nodes; //  32 * 2 - 1 
    bitstream_type*    bs;
} squeeze_type;

#define squeeze_size_mul(name, count) (                                \
    ((uint64_t)(count) >= ((SIZE_MAX / 4) / (uint64_t)sizeof(name))) ? \
    0 : (size_t)((uint64_t)sizeof(name) * (uint64_t)(count))           \
)

#define squeeze_size_implementation(win_bits) (                        \
    (sizeof(squeeze_type)) +                                           \
    /* lit_nodes: */                                                   \
    squeeze_size_mul(huffman_node_type, (512ULL * 2ULL - 1ULL)) +      \
    /* pos_nodes: */                                                   \
    squeeze_size_mul(huffman_node_type, ((1ULL << 5) * 2ULL - 1ULL))   \
)

#define squeeze_sizeof(win_bits) (                                     \
    (sizeof(size_t) == sizeof(uint64_t)) &&                            \
    (squeeze_min_win_bits <= (win_bits)) &&                            \
                             ((win_bits) <= squeeze_max_win_bits) ?    \
    (size_t)squeeze_size_implementation((win_bits)) : 0                \
)

typedef struct {
    errno_t (*init_with)(squeeze_type* s, void* memory, size_t size,
                         uint8_t win_bits);
    // `win_bits` is a log2 of window size in bytes in range
    // [squeeze_min_win_bits..squeeze_max_win_bits]
    void (*write_header)(bitstream_type* bs, uint64_t bytes, uint8_t win_bits);
    void (*compress)(squeeze_type* s, const uint8_t* data, size_t bytes, uint16_t window);
    void (*read_header)(bitstream_type* bs, uint64_t *bytes, uint8_t *win_bits);
    void (*decompress)(squeeze_type* s, uint8_t* data, size_t bytes);
} squeeze_interface;

extern squeeze_interface squeeze;

#if 0

// TODO: consider inclusion of these two functions for convenience:

static squeeze_type* squeeze_new(bitstream_type* bs, uint8_t win_bits,
                                 uint8_t map_bits, uint8_t len_bits) {
    const uint64_t bytes = squeeze_sizeof(win_bits, map_bits, len_bits);
    squeeze_type* s = (squeeze_type*)calloc(1, (size_t)bytes);
    if (s != null) {
        squeeze_init_with(s, s, bytes, win_bits, map_bits, len_bits);
        s->bs = bs;
    }
    return s;
}

static void squeeze_delete(squeeze_type* s) {
    free(s);
}

#endif

#endif // squeeze_header_included

#if defined(squeeze_implementation) && !defined(squeeze_implemented)

#define squeeze_implemented

#define bitstream_implementation
#include "bitstream.h"

#define huffman_implementation
#include "huffman.h"

#include "deflate.h"

#ifndef null
#define null ((void*)0) // like null_ptr better than NULL (0)
#endif

#ifndef countof
#define countof(a) (sizeof(a) / sizeof((a)[0]))
#endif

#ifndef assert
#include <assert.h>
#endif

static errno_t squeeze_init_with(squeeze_type* s, void* memory, size_t size,
                                 uint8_t win_bits) {
    errno_t r = 0;
    assert(squeeze_min_win_bits <= win_bits && win_bits <= squeeze_max_win_bits);
    size_t expected = squeeze_sizeof(win_bits);
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
        const size_t lit_m = lit_n * 2 - 1;
        const size_t pos_m = pos_n * 2 - 1;
        s->lit_nodes = (huffman_node_type*)p; p += sizeof(huffman_node_type) * lit_m;
        s->pos_nodes = (huffman_node_type*)p; p += sizeof(huffman_node_type) * pos_m;
        assert(p == (uint8_t*)memory + size);
        huffman_init(&s->lit, s->lit_nodes, lit_m);
        huffman_init(&s->pos, s->pos_nodes, pos_m);
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
        huffman_insert(&s->lit, lit);
    } else {
        squeeze_write_huffman(s, &s->lit, lit);
    }
}

static inline void squeeze_encode_len(squeeze_type* s, const uint16_t len) {
    assert(3 <= len && len <= 258);
    uint8_t  i = deflate_len_index[len];
    uint16_t b = deflate_len_base[i];
    uint8_t  x = deflate_len_extra_bits[i];
    squeeze_encode_literal(s, (uint16_t)(257 + i));
    assert(b <= len && len - b <= (uint16_t)(1u << x));
    if (x > 0) { squeeze_write_bits(s, (len - b), x); }
}

static inline void squeeze_encode_pos(squeeze_type* s, const uint16_t pos) {
    assert(0 <= pos && pos <= 0xFFFF);
    uint8_t  i = deflate_pos_index[pos];
    uint16_t b = deflate_pos_base[i];
    uint8_t  x = deflate_pos_extra_bits[i];
    if (s->pos.node[i].bits == 0) {
        assert(s->pos.node[squeeze_pos_nyt].bits != 0);
        squeeze_write_huffman(s, &s->pos, squeeze_pos_nyt);
        squeeze_write_bits(s, i, 5); // 0..29
        huffman_insert(&s->pos, i);
    } else {
        squeeze_write_huffman(s, &s->pos, i);
    }
    assert(b <= pos && pos - b <= (uint16_t)(1u << x));
    if (x > 0) { squeeze_write_bits(s, (pos - b), x); }
}

static void squeeze_compress(squeeze_type* s, const uint8_t* data, uint64_t bytes, 
                             uint16_t window) {
    size_t pl_bytes = 0;
    size_t as_bytes = 0;
    huffman_insert(&s->lit, squeeze_lit_nyt);
    huffman_insert(&s->pos, squeeze_pos_nyt);
    deflate_init();
    size_t i = 0;
    while (i < bytes && s->error == 0) {
        size_t len = 0;
        size_t pos = 0;
        if (i >= 1) {
            size_t j = i - 1;
            size_t min_j = i > window ? i - window : 0;
            while (j > min_j) {
                assert((i - j) < window);
                const size_t n = bytes - i;
                size_t k = 0;
                while (k < n && data[j + k] == data[i + k] && k < deflate_max_len) {
                    k++;
                }
                if (k > len) {
                    len = k;
                    pos = i - j;
                    if (len == deflate_max_len) { break; }
                }
                j--;
            }
        }
        if (len >= 3) {
            assert(0 < pos && pos < window);
            squeeze_encode_len(s, (uint16_t)len);
            squeeze_encode_pos(s, (uint16_t)pos);
            i += len;
            pl_bytes += len;
        } else {
            as_bytes++;
            squeeze_encode_literal(s, data[i]);
            i++;
        }
    }
    squeeze_flush(s);
#if 0
    double pl_percent = (100.0 * pl_bytes) / (pl_bytes + as_bytes);
    double as_percent = (100.0 * as_bytes) / (pl_bytes + as_bytes);
    printf("lit: %.2f%% back references: %.2f%% \n", as_percent, pl_percent);
    printf("entropy: lit: %.2f pos: %.2f\n",
        huffman_entropy(&s->lit), huffman_entropy(&s->pos));
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

static void squeeze_decompress(squeeze_type* s, uint8_t* data, uint64_t bytes) {
    huffman_insert(&s->lit, squeeze_lit_nyt);
    huffman_insert(&s->pos, squeeze_pos_nyt);
    deflate_init();
    size_t i = 0; // output b64[i]
    while (i < bytes && s->error == 0) {
        uint64_t lit = squeeze_read_huffman(s, &s->lit);
        if (lit == squeeze_lit_nyt) {
            lit = squeeze_read_bits(s, 9);
            huffman_insert(&s->lit, (int32_t)lit);
        }
        if (lit <= 0xFF) {
            data[i] = (uint8_t)lit;
            i++;
        } else {
            assert(257 <= lit && lit <= squeeze_lit_nyt);
            if (257 <= lit && lit <= squeeze_lit_nyt) {
                uint8_t  len_base   = (uint8_t)(lit - 257);
                uint8_t  len_bits   = deflate_len_extra_bits[len_base];
                uint64_t len_extra = len_bits == 0 ? 0 : squeeze_read_bits(s, len_bits);
                uint32_t len = (uint32_t)deflate_len_base[len_base] + (uint32_t)len_extra;
                uint64_t pos_base = squeeze_read_huffman(s, &s->pos);
                squeeze_if_error_return(s);
                if (pos_base == squeeze_pos_nyt) {
                    pos_base = squeeze_read_bits(s, 5);
                    squeeze_if_error_return(s);
                    huffman_insert(&s->pos, (int32_t)pos_base);
                }
                uint8_t  pos_bits  = deflate_pos_extra_bits[pos_base];
                uint64_t pos_extra = pos_bits == 0 ? 0 : squeeze_read_bits(s, pos_bits);
                uint32_t pos = (uint32_t)deflate_pos_base[pos_base] + (uint32_t)pos_extra;
                squeeze_if_error_return(s);
                assert(3 <= len);
                if (len < 3) { squeeze_return_invalid(s); }
                // Cannot do memcpy() here because of possible overlap.
                // memcpy() may read more than one byte at a time.
                uint8_t* d = data - (size_t)pos;
                const size_t n = i + (size_t)len;
                while (i < n) { data[i] = d[i]; i++; }
            } else {
                s->error = EINVAL;
            }
        }
    }
}

squeeze_interface squeeze = {
    .init_with    = squeeze_init_with,
    .write_header = squeeze_write_header,
    .compress     = squeeze_compress,
    .read_header  = squeeze_read_header,
    .decompress   = squeeze_decompress,
};

#endif // squeeze_implementation
