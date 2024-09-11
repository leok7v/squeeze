#ifndef squeeze_header_included
#define squeeze_header_included

#include <errno.h>
#include <stdint.h>

#include "bitstream.h"
#include "huffman.h"
#include "map.h"

enum {
    squeeze_syllable_bits = 4, // 4-bits syllables for "length" > 255
    squeeze_min_win_bits  = 10,
    squeeze_max_win_bits  = 20,
    squeeze_min_map_bits  =  8,
    squeeze_max_map_bits  = 20,
    squeeze_min_len_bits  =  4,
    squeeze_max_len_bits  =  8
};

typedef struct {
    errno_t error; // sticky
    map_type map;  // `words` dictionary
    map_entry_t* map_entries;
    huffman_tree_type dic; // `map` keys tree
    huffman_tree_type asc; // 0..127 ASCII characters
    huffman_tree_type bin; // 128..255 characters
    huffman_tree_type pos; // positions tree of 1^win_bits
    huffman_tree_type len; // length [2..255] tree
    huffman_node_type* dic_nodes;
    huffman_node_type* asc_nodes;
    huffman_node_type* bin_nodes;
    huffman_node_type* pos_nodes;
    huffman_node_type* len_nodes;
    bitstream_type*    bs;
} squeeze_type;

#define squeeze_size_mul(name, count) (                                         \
    ((uint64_t)(count) >= ((SIZE_MAX / 4) / (uint64_t)sizeof(name))) ?          \
    0 : (size_t)((uint64_t)sizeof(name) * (uint64_t)(count))                    \
)

#define squeeze_size_implementation(win_bits, map_bits, len_bits) (             \
    (sizeof(squeeze_type)) +                                                    \
    squeeze_size_mul(map_entry_t, (1ULL << (map_bits))) +                       \
    /* dic_nodes: */                                                            \
    squeeze_size_mul(huffman_node_type, ((1ULL << (map_bits)) * 2ULL - 1ULL)) + \
    /* asc_nodes: */                                                            \
    squeeze_size_mul(huffman_node_type, (128ULL * 2ULL - 1ULL)) +               \
    /* bin_nodes: */                                                            \
    squeeze_size_mul(huffman_node_type, (128ULL * 2ULL - 1ULL)) +               \
    /* pos_nodes: */                                                            \
    squeeze_size_mul(huffman_node_type, ((1ULL << (win_bits)) * 2ULL - 1ULL)) + \
    /* len_nodes: */                                                            \
    squeeze_size_mul(huffman_node_type, ((1ULL << (len_bits)) * 2ULL - 1ULL))   \
)

#define squeeze_sizeof(win_bits, map_bits, len_bits) (                          \
    (sizeof(size_t) == sizeof(uint64_t)) &&                                     \
    (squeeze_min_win_bits <= (win_bits)) &&                                     \
                             ((win_bits) <= squeeze_max_win_bits) &&            \
    (squeeze_min_map_bits <= (map_bits)) &&                                     \
                            ((map_bits) <= squeeze_max_map_bits) &&             \
    (squeeze_min_len_bits <= (len_bits)) &&                                     \
                            ((len_bits) <= squeeze_max_len_bits) ?              \
    (size_t)squeeze_size_implementation((win_bits), (map_bits), (len_bits)) : 0 \
)

typedef struct {
    errno_t (*init_with)(squeeze_type* s, void* memory, size_t size,
                         uint8_t win_bits, uint8_t map_bits,
                         uint8_t len_bits);
    // `win_bits` is a log2 of window size in bytes in range
    // [squeeze_min_win_bits..squeeze_max_win_bits]
    void (*write_header)(bitstream_type* bs, uint64_t bytes,
                         uint8_t win_bits, uint8_t map_bits, uint8_t len_bits);
    void (*compress)(squeeze_type* s, const uint8_t* data, size_t bytes);
    void (*read_header)(bitstream_type* bs, uint64_t *bytes,
                        uint8_t *win_bits, uint8_t *map_bits, uint8_t *len_bits);
    void (*decompress)(squeeze_type* s, uint8_t* data, size_t bytes);
} squeeze_interface;

extern squeeze_interface squeeze;

#if 0

// TODO: consider inclusion of these two functions for convinience:

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

#include "bitstream.h"

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
                                 uint8_t win_bits, uint8_t map_bits,
                                 uint8_t len_bits) {
    errno_t r = 0;
    assert(squeeze_min_win_bits <= win_bits && win_bits <= squeeze_max_win_bits);
    assert(squeeze_min_map_bits <= map_bits && map_bits <= squeeze_max_map_bits);
    assert(squeeze_min_len_bits <= len_bits && len_bits <= squeeze_max_len_bits);
    size_t expected = squeeze_sizeof(win_bits, map_bits, len_bits);
    // 167,936,192 bytes for (win_bits = 11, map_bits = 19)
    assert(size == expected);
    if (expected == 0 || memory == null || size != expected) {
        r = EINVAL;
    } else {
        uint8_t* p = (uint8_t*)memory;
        memset(memory, 0, sizeof(squeeze_type));
        p += sizeof(squeeze_type);
        const size_t map_n = ((size_t)1U) << map_bits;
        const size_t dic_n = map_n;
        const size_t asc_n = 128; // always 128
        const size_t bin_n = 128; // always 128
        const size_t pos_n = ((size_t)1U) << win_bits;
        const size_t len_n = ((size_t)1U) << len_bits;
        const size_t dic_m = dic_n * 2 - 1;
        const size_t asc_m = asc_n * 2 - 1;
        const size_t bin_m = bin_n * 2 - 1;
        const size_t pos_m = pos_n * 2 - 1;
        const size_t len_m = len_n * 2 - 1;
        s->map_entries = (map_entry_t*)p; p += sizeof(map_entry_t) * map_n;
        s->dic_nodes = (huffman_node_type*)p; p += sizeof(huffman_node_type) * dic_m;
        s->asc_nodes = (huffman_node_type*)p; p += sizeof(huffman_node_type) * asc_m;
        s->bin_nodes = (huffman_node_type*)p; p += sizeof(huffman_node_type) * bin_m;
        s->pos_nodes = (huffman_node_type*)p; p += sizeof(huffman_node_type) * pos_m;
        s->len_nodes = (huffman_node_type*)p; p += sizeof(huffman_node_type) * len_m;
        assert(p == (uint8_t*)memory + size);
        map.init(&s->map,     s->map_entries, map_n);
        huffman.init(&s->asc, s->asc_nodes, asc_m);
        huffman.init(&s->bin, s->bin_nodes, bin_m);
        huffman.init(&s->dic, s->dic_nodes, dic_m);
        huffman.init(&s->pos, s->pos_nodes, pos_m);
        huffman.init(&s->len, s->len_nodes, len_m);
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
        bitstream.write_bit(s->bs, bit);
        s->error = s->bs->error;
    }
}

static inline void squeeze_write_bits(squeeze_type* s,
                                      uint64_t b64, uint8_t bits) {
    if (s->error == 0) {
        bitstream.write_bits(s->bs, b64, bits);
        s->error = s->bs->error;
    }
}

static inline int32_t squeeze_bit_count(uint64_t number, uint8_t base) {
    assert(number > 0);
    int32_t bits = 0;
    while (number != 0) {
        bits += base;
        number >>= base;
        bits++; // continue bit
    }
    return bits;
}

static inline int32_t squeeze_len_pos_bits(const squeeze_type* s,
                                           size_t len, size_t pos,
                                           uint8_t len_bits) {
    if (len < 3) {
        return 0;
    } if (len < (1ULL << len_bits)) {
        return s->len.node[len].bits + s->pos.node[pos].bits;
    } else {
        return s->len.node[0].bits +
               squeeze_bit_count(len, squeeze_syllable_bits) +
               s->pos.node[pos].bits;
    }
}

static inline void squeeze_write_number(squeeze_type* s,
                                        uint64_t bits, uint8_t base) {
    while (s->error == 0 && bits != 0) {
        squeeze_write_bits(s, bits, base);
        bits >>= base;
        squeeze_write_bit(s, bits != 0); // continue bit
    }
}

static inline void squeeze_write_huffman(squeeze_type* s, huffman_tree_type* t,
                                         int32_t i) {
    assert(t != null && t->node != null);
    assert(0 <= i && i < t->n); // leaf symbol
    assert(1 <= t->node[i].bits && t->node[i].bits < 64);
    squeeze_write_bits(s, t->node[i].path, (uint8_t)t->node[i].bits);
    huffman.inc_frequency(t, i); // after the path is written
}

static inline void squeeze_flush(squeeze_type* s) {
    if (s->error == 0) {
        bitstream.flush(s->bs);
        s->error = s->bs->error;
    }
}

static void squeeze_write_header(bitstream_type* bs, uint64_t bytes,
                                 uint8_t win_bits, uint8_t map_bits,
                                 uint8_t len_bits) {
    if (win_bits < squeeze_min_win_bits || win_bits > squeeze_max_win_bits ||
        map_bits < squeeze_min_map_bits || map_bits > squeeze_max_map_bits ||
        len_bits < squeeze_min_len_bits || len_bits > squeeze_max_len_bits) {
        bs->error = EINVAL;
    } else {
        enum { bits64 = sizeof(uint64_t) * 8 };
        bitstream.write_bits(bs, (uint64_t)bytes, bits64);
        enum { bits8 = sizeof(uint8_t) * 8 };
        bitstream.write_bits(bs, win_bits, bits8);
        bitstream.write_bits(bs, map_bits, bits8);
        bitstream.write_bits(bs, len_bits, bits8);
    }
}

static void squeeze_add_to_dictionary(squeeze_type* s, const uint8_t* word,
                                      uint64_t bytes) {
    size_t word_bytes = bytes < countof(s->map_entries[0]) - 1 ?
                        bytes : countof(s->map_entries[0]) - 1;
    assert(word_bytes <= 0xFF);
    int32_t wix = map.put(&s->map, word, (uint8_t)word_bytes);
    if (wix >= 0) {
        huffman.inc_frequency(&s->dic, wix);
    }
}

static void squeeze_compress(squeeze_type* s, const uint8_t* data, uint64_t bytes) {
size_t dc_bytes = 0;
size_t pl_bytes = 0;
size_t as_bytes = 0;
    squeeze_if_error_return(s);
    const uint8_t win_bits = huffman.log2_of_pow2(s->pos.n);
    const uint8_t len_bits = huffman.log2_of_pow2(s->len.n);
    if (win_bits < 10 || win_bits > 20) { squeeze_return_invalid(s); }
    const size_t window = ((size_t)1U) << win_bits;
    size_t i = 0;
    while (i < bytes) {
        // "best" longest dictionary match
        const int32_t best = map.best(&s->map, &data[i], bytes - i);
        const uint8_t best_bytes = best < 0 ? 0 : map.bytes(&s->map, best);
        const int32_t best_bits  = best < 0 ? 0 :
                        s->dic.node[best].bits + s->len.node[1].bits;
        // bytes and position of longest matching sequence
        size_t len = 0;
        size_t pos = 0;
        if (i >= 1) {
            size_t j = i - 1;
            size_t min_j = i > window ? i - window : 0;
            while (j > min_j) {
                assert((i - j) < window);
                const size_t n = bytes - i;
                size_t k = 0;
                while (k < n && data[j + k] == data[i + k]) {
                    k++;
                }
                if (k > len) {
                    len = k;
                    pos = i - j;
                }
                j--;
            }
        }
        const uint32_t len_pos_bits = squeeze_len_pos_bits(s, len, pos,
                                                           len_bits);
        if (0 < len_pos_bits && len_pos_bits < len * 8) {
            assert(0 < pos && pos < window);
            squeeze_write_bits(s, 0b11, 2); // flags
            squeeze_if_error_return(s);
            if (len < (1ULL << len_bits)) {
                squeeze_write_huffman(s, &s->len, (int32_t)len);
            } else {
                squeeze_write_huffman(s, &s->len, 0);
                squeeze_if_error_return(s);
                squeeze_write_number(s, len, squeeze_syllable_bits);
            }
            squeeze_if_error_return(s);
            squeeze_write_huffman(s, &s->pos, (int32_t)pos);
            squeeze_if_error_return(s);
            squeeze_add_to_dictionary(s, &data[i], len);
            i += len;
            pl_bytes += len;
        } else {
            if (best_bits < best_bytes * 8) {
                assert(best <= INT32_MAX);
                assert(map.bytes(&s->map, best) >= 3);
                squeeze_write_bits(s, 0b11, 2); // flags
                squeeze_if_error_return(s);
                // len == 1 indicates that it's a dictionary word
                squeeze_write_huffman(s, &s->len, 1);
                squeeze_if_error_return(s);
                assert(s->dic.node[best].bits <= 0xFF);
                squeeze_write_huffman(s, &s->dic, (int32_t)best);
                squeeze_if_error_return(s);
                i += best_bytes;
                dc_bytes += best_bytes;
            } else {
                as_bytes++;
                const uint8_t b = data[i];
                // European texts are predominantly spaces and small ASCII letters:
                if (b < 0x80) {
                    squeeze_write_bit(s, 0); // flags
                    squeeze_if_error_return(s);
                    // ASCII byte < 0x80 with 8th bit set to `0`
                    squeeze_write_huffman(s, &s->asc, b);
                    squeeze_if_error_return(s);
                } else {
                    squeeze_write_bit(s, 1); // flag: 1
                    squeeze_write_bit(s, 0); // flag: 0
                    squeeze_if_error_return(s);
                    // only 7 bit because 8th bit is `1`
                    squeeze_write_huffman(s, &s->bin, b & 0x7F);
                    squeeze_if_error_return(s);
                }
                i++;
            }
        }
    }
    squeeze_flush(s);
//  printf("dic: %lld p/l: %lld sym: %lld\n", dc_bytes, pl_bytes, as_bytes);
    double dc_percent = (100.0 * dc_bytes) / (dc_bytes + pl_bytes + as_bytes);
    double pl_percent = (100.0 * pl_bytes) / (dc_bytes + pl_bytes + as_bytes);
    double as_percent = (100.0 * as_bytes) / (dc_bytes + pl_bytes + as_bytes);
    printf("dic: %.2f%% p/l: %.2f%% sym: %.2f%%\n", dc_percent, pl_percent, as_percent);
    printf("entropy: dic: %.2f len: %.2f pos: %.2f asc: %.2f bin: %.2f\n",
        huffman.entropy(&s->dic), huffman.entropy(&s->len),
        huffman.entropy(&s->pos), huffman.entropy(&s->asc),
        huffman.entropy(&s->bin));
//   printf("s->len[0].bits: %d s->len[1].bits: %d\n", s->len.node[0].bits, s->len.node[1].bits);
}

static inline uint64_t squeeze_read_bit(squeeze_type* s) {
    bool bit = 0;
    if (s->error == 0) {
        bit = bitstream.read_bit(s->bs);
        s->error = s->bs->error;
    }
    return bit;
}

static inline uint64_t squeeze_read_bits(squeeze_type* s, uint32_t n) {
    assert(n <= 64);
    uint64_t bits = 0;
    if (s->error == 0) {
        bits = bitstream.read_bits(s->bs, n);
    }
    return bits;
}

static inline uint64_t squeeze_read_number(squeeze_type* s, uint8_t base) {
    uint64_t bits = 0;
    uint64_t bit = 0;
    uint32_t shift = 0;
    while (s->error == 0) {
        bits |= (squeeze_read_bits(s, base) << shift);
        shift += base;
        bit = squeeze_read_bit(s);
        if (!bit) { break; }
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
    assert(0 <= i && i < t->n); // leaf symbol
    huffman.inc_frequency(t, i);
    return (uint64_t)i;
}

static void squeeze_read_header(bitstream_type* bs, uint64_t *bytes,
                                uint8_t *win_bits, uint8_t *map_bits,
                                uint8_t *len_bits) {
    uint64_t b  = bitstream.read_bits(bs, sizeof(uint64_t) * 8);
    uint64_t wb = bitstream.read_bits(bs, sizeof(uint8_t) * 8);
    uint64_t mb = bitstream.read_bits(bs, sizeof(uint8_t) * 8);
    uint64_t lb = bitstream.read_bits(bs, sizeof(uint8_t) * 8);
    if (bs->error == 0) {
        if (wb < squeeze_min_win_bits || wb > squeeze_max_win_bits) {
            bs->error = EINVAL;
        } else if (mb < squeeze_min_map_bits || mb > squeeze_max_map_bits) {
            bs->error = EINVAL;
        } else if (lb < 4 || lb > 8) {
            bs->error = EINVAL;
        } else if (bs->error == 0) {
            *bytes = b;
            *win_bits = (uint8_t)wb;
            *map_bits = (uint8_t)mb;
            *len_bits = (uint8_t)lb;
        }
    }
}

static void squeeze_decompress(squeeze_type* s, uint8_t* data, uint64_t bytes) {
    squeeze_if_error_return(s);
    const uint8_t win_bits = huffman.log2_of_pow2(s->pos.n);
    if (win_bits < 10 || win_bits > 20) { squeeze_return_invalid(s); }
    const size_t window = ((size_t)1U) << win_bits;
    size_t i = 0; // output b64[i]
    while (i < bytes) {
        uint64_t bit0 = squeeze_read_bit(s);
        squeeze_if_error_return(s);
        if (bit0) {
            uint64_t bit1 = squeeze_read_bit(s);
            squeeze_if_error_return(s);
            if (bit1) {
                uint64_t len = squeeze_read_huffman(s, &s->len);
                squeeze_if_error_return(s);
                if (len == 1) {
                    uint64_t wix = squeeze_read_huffman(s, &s->dic);
                    squeeze_if_error_return(s);
                    assert(wix < (uint64_t)s->map.n);
                    assert(map.bytes(&s->map, (int32_t)wix) >= 3);
                    size_t n = map.bytes(&s->map, (int32_t)wix);
                    assert(i + n <= bytes);
                    const uint8_t* d = (const uint8_t*)map.data(&s->map, (int32_t)wix);
                    for (size_t j = 0; j < n; j++) { data[i] = d[j]; i++; }
                } else {
                    if (len == 0) { len = squeeze_read_number(s, squeeze_syllable_bits); }
                    uint64_t pos = squeeze_read_huffman(s, &s->pos);
                    squeeze_if_error_return(s);
                    assert(0 < pos && pos < window);
                    if (!(0 < pos && pos < window)) { squeeze_return_invalid(s); }
                    assert(2 <= len);
                    if (len < 2) { squeeze_return_invalid(s); }
                    // Cannot do memcpy() here because of possible overlap.
                    // memcpy() may read more than one byte at a time.
                    uint8_t* d = data - (size_t)pos;
                    const size_t n = i + (size_t)len;
                    uint8_t* w = d + i;
                    while (i < n) { data[i] = d[i]; i++; }
                    squeeze_add_to_dictionary(s, w, len);
                }
            } else { // byte >= 0x80
                uint64_t b = squeeze_read_huffman(s, &s->bin);
                squeeze_if_error_return(s);
                data[i] = (uint8_t)b | 0x80;
                i++;
            }
        } else { // literal byte (ASCII byte < 0x80)
            uint64_t b = squeeze_read_huffman(s, &s->asc);
            squeeze_if_error_return(s);
            data[i] = (uint8_t)b;
            i++;
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
