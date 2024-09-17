#ifndef squeeze_header_included
#define squeeze_header_included

#include <errno.h>
#include <stdint.h>

enum {
    squeeze_deflate_sym_min   = 257, // minimum literal for length base
    squeeze_deflate_sym_max   = 284, // maximum literal for length base
    squeeze_deflate_pos_max   = 29,   // maximum pos base index
    squeeze_deflate_len_min   = 3,
    // value is the same but unrelated to squeeze_deflate_sym_min:
    squeeze_deflate_len_max   = 257
};

enum { // "nyt" stands for Not Yet Transmitted (see Vitter Algorithm)
    squeeze_min_win_bits  =  10,
    squeeze_max_win_bits  =  15,
    squeeze_lit_nyt       =  squeeze_deflate_sym_max + 1,
    squeeze_pos_nyt       =  squeeze_deflate_pos_max + 1,
};

// deflate tables

static const uint16_t squeeze_len_base[29] = {
    3, 4, 5, 6, 7, 8, 9, 10, // 257-264
    11, 13, 15, 17,          // 265-268
    19, 23, 27, 31,          // 269-272
    35, 43, 51, 59,          // 273-276
    67, 83, 99, 115,         // 277-280
    131, 163, 195, 227, 258  // 281-285
};

static const uint8_t squeeze_len_xb[29] = { // extra bits
    0, 0, 0, 0, 0, 0, 0, 0, // 257-264
    1, 1, 1, 1,             // 265-268
    2, 2, 2, 2,             // 269-272
    3, 3, 3, 3,             // 273-276
    4, 4, 4, 4,             // 277-280
    5, 5, 5, 5, 0           // 281-285 (len = 258 has no extra bits)
};

static const uint16_t squeeze_pos_base[30] = {
    1, 2, 3, 4,  // 0-3
    5, 7,        // 4-5
    9, 13,       // 6-7
    17, 25,      // 8-9
    33, 49,      // 10-11
    65, 97,      // 12-13
    129, 193,    // 14-15
    257, 385,    // 16-17
    513, 769,    // 18-19
    1025, 1537,  // 20-21
    2049, 3073,  // 22-23
    4097, 6145,  // 24-25
    8193, 12289, // 26-27
    16385, 24577 // 28-29
};

static const uint8_t squeeze_pos_xb[30] = { // extra bits
    0, 0, 0, 0, // 0-3
    1, 1,       // 4-5
    2, 2,       // 6-7
    3, 3,       // 8-9
    4, 4,       // 10-11
    5, 5,       // 12-13
    6, 6,       // 14-15
    7, 7,       // 16-17
    8, 8,       // 18-19
    9, 9,       // 20-21
    10, 10,     // 22-23
    11, 11,     // 24-25
    12, 12,     // 26-27
    13, 13      // 28-29
};

typedef struct huffman_node {
    uint64_t freq;
    uint64_t path;
    int32_t  bits; // 0 for root
    int32_t  pix;  // parent
    int32_t  lix;  // left
    int32_t  rix;  // right
} huffman_node;

typedef struct huffman_tree {
    huffman_node* node;
    int32_t n;
    int32_t next;  // next non-terminal nodes in the tree >= n
    int32_t depth; // max tree depth seen
    int32_t complete; // tree is too deep or freq too high - no more updates
    // stats:
    struct {
        size_t updates;
        size_t swaps;
        size_t moves;
    } stats;
} huffman_tree;

typedef struct bitstream {
    void*    stream; // stream and (data,capacity) are exclusive
    uint8_t* data;
    uint64_t capacity; // data[capacity]
    uint64_t bytes; // number of bytes written
    uint64_t read;  // number of bytes read
    uint64_t b64;   // bit shifting buffer
    int32_t  bits;  // bit count inside b64
    errno_t  error; // sticky error
    errno_t (*output)(struct bitstream* bs); // write b64 as 8 bytes
    errno_t (*input)(struct bitstream* bs);  // read b64 as 8 bytes
} bitstream;

typedef struct {
    errno_t error; // sticky
    huffman_tree lit; // 0..255 literal bytes; 257-285 length
    huffman_tree pos; // positions tree of 1^win_bits
    huffman_node lit_nodes[512 * 2 - 1];
    huffman_node pos_nodes[32 * 2 - 1];
    bitstream*         bs;
    uint8_t len_index[squeeze_deflate_sym_max + 1]; // index squeeze_len_base
    uint8_t pos_index[1u << 15];                    // index squeeze_pos_base
} squeeze_type;

typedef struct {
    void (*init)(squeeze_type* s);
    // `win_bits` is a log2 of window size in bytes in range
    // [squeeze_min_win_bits..squeeze_max_win_bits]
    void (*write_header)(bitstream* bs, uint64_t bytes);
    void (*compress)(squeeze_type* s, bitstream* bs,
                     const uint8_t* data, size_t bytes, uint16_t window);
    void (*read_header)(bitstream* bs, uint64_t *bytes);
    void (*decompress)(squeeze_type* s, bitstream* bs,
                       uint8_t* data, size_t bytes);
} squeeze_interface;

#if defined(__cplusplus)
extern "C" squeeze_interface squeeze;
#else
extern squeeze_interface squeeze;
#endif

#endif // squeeze_header_included

#if defined(squeeze_implementation) && !defined(squeeze_implemented)

#define squeeze_implemented

// https://stackoverflow.com/questions/76899757/how-can-i-use-nullptr-for-null-pointers-prior-to-c23

#ifndef null
#define null ((void*)0) // like null_ptr better than NULL (0)
#endif

#ifndef countof
#define countof(a) (sizeof(a) / sizeof((a)[0]))
#endif

#ifndef assert
#include <assert.h>
#endif
#include <string.h>

static inline void     bitstream_create(bitstream* bs, void* data,
                                                       size_t capacity);
static inline void     bitstream_write_bit(bitstream* bs, int32_t bit);
static inline void     bitstream_write_bits(bitstream* bs, uint64_t data,
                                                           int32_t bits);
static inline int      bitstream_read_bit(bitstream* bs); // bool 0|1 "int"
static inline uint64_t bitstream_read_bits(bitstream* bs, int32_t bits);
static inline void     bitstream_flush(bitstream* bs); // write trailing zeros
static inline void     bitstream_dispose(bitstream* bs);

static inline void bitstream_write_bit(bitstream* bs, int32_t bit) {
    if (bs->error == 0) {
        bs->b64 <<= 1;
        bs->b64 |= (bit & 1);
        bs->bits++;
        if (bs->bits == 64) {
            if (bs->data != null && bs->capacity > 0) {
                assert(bs->stream == null);
                for (int i = 0; i < 8 && bs->error == 0; i++) {
                    if (bs->bytes == bs->capacity) {
                        bs->error = E2BIG;
                    } else {
                        uint8_t b = (uint8_t)(bs->b64 >> ((7 - i) * 8));
                        bs->data[bs->bytes++] = b;
                    }
                }
            } else {
                assert(bs->data == null && bs->capacity == 0);
                bs->error = bs->output(bs);
                if (bs->error == 0) { bs->bytes += 8; }
            }
            bs->bits = 0;
            bs->b64 = 0;
        }
    }
}

static inline void bitstream_write_bits(bitstream* bs,
                                        uint64_t data, int32_t bits) {
    assert(0 < bits && bits <= 64);
    while (bits > 0 && bs->error == 0) {
        bitstream_write_bit(bs, data & 1);
        bits--;
        data >>= 1;
    }
}

static inline int bitstream_read_bit(bitstream* bs) {
    int bit = 0;
    if (bs->error == 0) {
        if (bs->bits == 0) {
            bs->b64 = 0;
            if (bs->data != null && bs->bytes > 0) {
                assert(bs->stream == null);
                for (int i = 0; i < 8 && bs->error == 0; i++) {
                    if (bs->read == bs->bytes) {
                        bs->error = E2BIG;
                    } else {
                        const uint64_t byte = (bs->data[bs->read] & 0xFF);
                        bs->b64 |= byte << ((7 - i) * 8);
                        bs->read++;
                    }
                }
            } else {
                assert(bs->data == null && bs->bytes == 0);
                bs->error = bs->input(bs);
                if (bs->error == 0) { bs->read += 8; }
            }
            bs->bits = 64;
        }
        bit = ((int64_t)bs->b64) < 0; // same as (bs->b64 >> 63) & 1;
        bs->b64 <<= 1;
        bs->bits--;
    }
    return bit;
}

static inline uint64_t bitstream_read_bits(bitstream* bs, int32_t bits) {
    uint64_t data = 0;
    assert(0 < bits && bits <= 64);
    for (int32_t b = 0; b < bits && bs->error == 0; b++) {
        int bit = bitstream_read_bit(bs);
        if (bit) { data |= ((uint64_t)bit) << b; }
    }
    return data;
}

static inline void bitstream_create(bitstream* bs, void* data, size_t bytes) {
    assert(bs->data != null);
    memset(bs, 0x00, sizeof(*bs));
    bs->data = (uint8_t*)data;
    bs->capacity = bytes;
}

static inline void bitstream_flush(bitstream* bs) {
    while (bs->bits > 0 && bs->error == 0) { bitstream_write_bit(bs, 0); }
}

static inline void bitstream_dispose(bitstream* bs) {
    memset(bs, 0x00, sizeof(*bs));
}

// Huffman Adaptive Coding https://en.wikipedia.org/wiki/Adaptive_Huffman_coding

static inline void   huffman_init(huffman_tree* t, huffman_node nodes[],
                                                   const size_t m);
static inline bool   huffman_insert(huffman_tree* t, int32_t i);
static inline bool   huffman_inc_frequency(huffman_tree* t, int32_t symbol);
static inline double huffman_entropy(const huffman_tree* t); // Shannon (bps)

static void huffman_update_paths(huffman_tree* t, int32_t i) {
    t->stats.updates++;
    const int32_t m = t->n * 2 - 1;
    if (i == m - 1) { t->depth = 0; } // root
    const int32_t  bits = t->node[i].bits;
    const uint64_t path = t->node[i].path;
    assert(bits < (int32_t)sizeof(uint64_t) * 8 - 1);
    assert((path & (~((1ULL << (bits + 1)) - 1))) == 0);
    const int32_t lix = t->node[i].lix;
    const int32_t rix = t->node[i].rix;
    if (lix != -1) {
        t->node[lix].bits = bits + 1;
        t->node[lix].path = path;
        huffman_update_paths(t, lix);
    }
    if (rix != -1) {
        t->node[rix].bits = bits + 1;
        t->node[rix].path = path | (1ULL << bits);
        huffman_update_paths(t, rix);
    }
    if (bits > t->depth) { t->depth = bits; }
}

static inline int32_t huffman_swap_siblings(huffman_tree* t,
                                            const int32_t i) {
    const int32_t m = t->n * 2 - 1;
    assert(0 <= i && i < m);
    if (i < m - 1) { // not root
        const int32_t pix = t->node[i].pix;
        assert(pix >= t->n); // parent (cannot be a leaf)
        const int32_t lix = t->node[pix].lix;
        const int32_t rix = t->node[pix].rix;
        if (lix >= 0 && rix >= 0) {
            assert(0 <= lix && lix < m - 1 && 0 <= rix && rix < m - 1);
            if (t->node[lix].freq > t->node[rix].freq) { // swap
                t->stats.swaps++;
                t->node[pix].lix = rix;
                t->node[pix].rix = lix;
                // because swap changed all path below:
                huffman_update_paths(t, pix);
                return i == lix ? rix : lix;
            }
        }
    }
    return i;
}

static inline void huffman_frequency_changed(huffman_tree* t, int32_t ix);

static inline void huffman_update_freq(huffman_tree* t, int32_t i) {
    const int32_t lix = t->node[i].lix;
    const int32_t rix = t->node[i].rix;
    assert(lix != -1 || rix != -1); // at least one leaf present
    t->node[i].freq = (lix >= 0 ? t->node[lix].freq : 0) +
                      (rix >= 0 ? t->node[rix].freq : 0);
}

static inline void huffman_move_up(huffman_tree* t, int32_t ix) {
    const int32_t pix = t->node[ix].pix; // parent
    assert(pix != -1);
    const int32_t gix = t->node[pix].pix; // grandparent
    assert(gix != -1);
    assert(t->node[pix].rix == ix);
    // Is parent grandparent`s left or right child?
    const bool parent_is_left_child = pix == t->node[gix].lix;
    const int32_t psx = parent_is_left_child ? // parent sibling index
        t->node[gix].rix : t->node[gix].lix;   // aka auntie/uncle
    if (t->node[ix].freq > t->node[psx].freq) {
        // Move grandparents left or right subtree to be
        // parents right child instead of 'i'.
        t->stats.moves++;
        t->node[ix].pix = gix;
        if (parent_is_left_child) {
            t->node[gix].rix = ix;
        } else {
            t->node[gix].lix = ix;
        }
        t->node[pix].rix = psx;
        t->node[psx].pix = pix;
        huffman_update_freq(t, pix);
        huffman_update_freq(t, gix);
        huffman_swap_siblings(t, ix);
        huffman_swap_siblings(t, psx);
        huffman_swap_siblings(t, pix);
        huffman_update_paths(t, gix);
        huffman_frequency_changed(t, gix);
    }
}

static inline void huffman_frequency_changed(huffman_tree* t, int32_t i) {
    const int32_t m = t->n * 2 - 1; (void)m;
    const int32_t pix = t->node[i].pix;
    if (pix == -1) { // `i` is root
        assert(i == m - 1);
        huffman_update_freq(t, i);
        i = huffman_swap_siblings(t, i);
    } else {
        assert(0 <= pix && pix < m);
        huffman_update_freq(t, pix);
        i = huffman_swap_siblings(t, i);
        huffman_frequency_changed(t, pix);
    }
    if (pix != -1 && t->node[pix].pix != -1 && i == t->node[pix].rix) {
        assert(t->node[i].freq >= t->node[t->node[pix].lix].freq);
        huffman_move_up(t, i);
    }
}

static bool huffman_insert(huffman_tree* t, int32_t i) {
    bool done = true;
    const int32_t root = t->n * 2 - 1 - 1;
    int32_t ipx = root;
    assert(t->node[i].pix == -1 && t->node[i].lix == -1 && t->node[i].rix == -1);
    assert(t->node[i].freq == 0 && t->node[i].bits == 0 && t->node[i].path == 0);
    t->node[i].freq = 1;
    while (ipx >= t->n) {
        if (t->node[ipx].rix == -1) {
            t->node[ipx].rix = i;
            t->node[i].pix = ipx;
            break;
        } else if (t->node[ipx].lix == -1) {
            t->node[ipx].lix = i;
            t->node[i].pix = ipx;
            break;
        } else {
            assert(t->node[ipx].lix >= 0);
            assert(t->node[i].freq <= t->node[t->node[ipx].lix].freq);
            ipx = t->node[ipx].lix;
        }
    }
    if (ipx >= t->n) { // not a leaf, inserted
        t->node[ipx].freq++;
        i = huffman_swap_siblings(t, i);
        assert(t->node[ipx].lix == i || t->node[ipx].rix);
        assert(t->node[ipx].freq ==
                (t->node[ipx].rix >= 0 ? t->node[t->node[ipx].rix].freq : 0) +
                (t->node[ipx].lix >= 0 ? t->node[t->node[ipx].lix].freq : 0));
    } else { // leaf
        assert(t->next > t->n);
        if (t->next == t->n) {
            done = false; // cannot insert
            t->complete = true;
        } else {
            t->next--;
            int32_t nix = t->next;
            t->node[nix] = (huffman_node){
                .freq = t->node[ipx].freq,
                .lix = ipx,
                .rix = -1,
                .pix = t->node[ipx].pix,
                .bits = t->node[ipx].bits,
                .path = t->node[ipx].path
            };
            if (t->node[ipx].pix != -1) {
                if (t->node[t->node[ipx].pix].lix == ipx) {
                    t->node[t->node[ipx].pix].lix = nix;
                } else {
                    t->node[t->node[ipx].pix].rix = nix;
                }
            }
            t->node[ipx].pix = nix;
            t->node[ipx].bits++;
            t->node[ipx].path = t->node[nix].path;
            t->node[nix].rix = i;
            t->node[i].pix = nix;
            t->node[i].bits = t->node[nix].bits + 1;
            t->node[i].path = t->node[nix].path | (1ULL << t->node[nix].bits);
            huffman_update_freq(t, nix);
            ipx = nix;
        }
    }
    huffman_frequency_changed(t, i);
    huffman_update_paths(t, ipx);
    assert(t->node[i].freq != 0 && t->node[i].bits != 0);
    return done;
}

// If input literal frequencies are severely skewed (e.g. Lucas numbers
// similar to Fibonacci numbers) and input sequence is long enough.
// The depth of the tree will grow past 64 bits.
// The first Lucas number that exceeds 2^64 is
// L(81) = 18,446,744,073,709,551,616 not actually realistic but
// better be safe than sorry

// huffman_inc_frequency() returns false on conditions inability to insert
// new node and tree being too deep of frequency too high
// The inability to insert new node should be treated as fatal outside
// huffman_inc_frequency() call.

static inline bool huffman_inc_frequency(huffman_tree* t, int32_t i) {
    bool done = true;
    assert(0 <= i && i < t->n); // terminal
    if (t->node[i].pix == -1) {
        done = huffman_insert(t, i); // Unseen terminal node.
    } else if (!t->complete && t->depth < 63 && t->node[i].freq < UINT64_MAX - 1) {
        t->node[i].freq++;
        huffman_frequency_changed(t, i);
    } else {
        // ignore future frequency updates
        t->complete = 1;
        done = false;
    }
    return done;
}

static inline double huffman_entropy(const huffman_tree* t) {
    // Shannon entropy
    double total = 0;
    double entropy = 0.0;
    for (int32_t i = 0; i < t->n; i++) { total += (double)t->node[i].freq; }
    for (int32_t i = 0; i < t->n; i++) {
        if (t->node[i].freq > 0) {
            double p_i = (double)t->node[i].freq / total;
            entropy += p_i * log2(p_i);
        }
    }
    return -entropy;
}

static inline void huffman_init(huffman_tree* t,
                                huffman_node nodes[],
                                const size_t count) {
    assert(7 <= count && count < INT32_MAX); // must pow(2, bits_per_symbol) * 2 - 1
    const int32_t n = (int32_t)(count + 1) / 2;
    assert(n > 4 && (n & (n - 1)) == 0); // must be power of 2
    memset(&t->stats, 0x00, sizeof(t->stats));
    t->node = nodes;
    t->n = n;
    const int32_t root = n * 2 - 1;
    t->next = root - 1; // next non-terminal node
    t->depth = 0;
    t->complete = 0;
    for (size_t i = 0; i < count; i++) {
        t->node[i] = (huffman_node){
            .freq = 0, .pix = -1, .lix = -1, .rix = -1, .bits = 0, .path = 0
        };
    }
}

// Deflate https://en.wikipedia.org/wiki/Deflate

static void squeeze_deflate_init(squeeze_type* s) {
    uint8_t  j = 0;
    uint16_t n = squeeze_len_base[j] + (1u << squeeze_len_xb[j]);
    for (int16_t i = 3; i < countof(s->len_index); i++) {
        if (i == n) {
            j++;
            n = squeeze_len_base[j] + (1u << squeeze_len_xb[j]);
        }
        s->len_index[i] = j;
    }
    assert(j == countof(squeeze_len_base) - 1);
    j = 0;
    n = squeeze_pos_base[j] + (1u << squeeze_pos_xb[j]);
    for (int32_t i = 0; i < countof(s->pos_index); i++) {
        if (i == n) {
            j++;
            n = squeeze_pos_base[j] + (1u << squeeze_pos_xb[j]);
        }
        s->pos_index[i] = j;
    }
    assert(j == countof(squeeze_pos_base) - 1);
}

// Squeeze https://en.wikipedia.org/wiki/LZ77_and_LZ78

static void squeeze_init(squeeze_type* s) {
    huffman_init(&s->lit, s->lit_nodes, countof(s->lit_nodes));
    huffman_init(&s->pos, s->pos_nodes, countof(s->pos_nodes));
}

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

static inline void squeeze_write_huffman(squeeze_type* s, huffman_tree* t,
                                         int32_t i) {
    assert(t != null && t->node != null);
    assert(0 <= i && i < t->n); // leaf symbol (literal)
    assert(1 <= t->node[i].bits && t->node[i].bits < 64);
    squeeze_write_bits(s, t->node[i].path, (uint8_t)t->node[i].bits);
    (void)huffman_inc_frequency(t, i); // after the path is written
    // (void) because inability to update frequency is OK here
}

static inline void squeeze_flush(squeeze_type* s) {
    if (s->error == 0) {
        bitstream_flush(s->bs);
        s->error = s->bs->error;
    }
}

static void squeeze_write_header(bitstream* bs, uint64_t bytes) {
    enum { bits64 = sizeof(uint64_t) * 8 };
    bitstream_write_bits(bs, (uint64_t)bytes, bits64);
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
    assert(0 <= lit && lit <= squeeze_deflate_sym_max);
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
    assert(3 <= len && len < sizeof(s->len_index));
    uint8_t  i = s->len_index[len];
    uint16_t b = squeeze_len_base[i];
    uint8_t  x = squeeze_len_xb[i];
    squeeze_encode_literal(s, (uint16_t)(squeeze_deflate_sym_min + i));
    assert(b <= len && len - b <= (uint16_t)(1u << x));
    if (x > 0) { squeeze_write_bits(s, (len - b), x); }
}

static inline void squeeze_encode_pos(squeeze_type* s, const uint16_t pos) {
    assert(0 <= pos && pos <= 0x7FFF);
    uint8_t  i = s->pos_index[pos];
    uint16_t b = squeeze_pos_base[i];
    uint8_t  x = squeeze_pos_xb[i];
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

static void squeeze_compress(squeeze_type* s, bitstream* bs,
                             const uint8_t* data, uint64_t bytes,
                             uint16_t window) {
    s->bs = bs;
    if (!huffman_insert(&s->lit, squeeze_lit_nyt)) { s->error = EINVAL; }
    if (!huffman_insert(&s->pos, squeeze_pos_nyt)) { s->error = EINVAL; }
    squeeze_deflate_init(s);
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
                while (k < n && data[j + k] == data[i + k] && k < squeeze_deflate_len_max) {
                    k++;
                }
                if (k >= squeeze_deflate_len_min && k > len) {
                    len = k;
                    pos = i - j;
                    if (len == squeeze_deflate_len_max) { break; }
                }
                if (j == min_j) { break; }
                j--;
            }
        }
        if (len >= squeeze_deflate_len_min) {
            assert(0 < pos && pos <= 0x7FFF);
            squeeze_encode_len(s, (uint16_t)len);
            squeeze_encode_pos(s, (uint16_t)pos);
            i += len;
        } else {
            squeeze_encode_literal(s, data[i]);
            i++;
        }
    }
    squeeze_flush(s);
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

static inline uint64_t squeeze_read_huffman(squeeze_type* s, huffman_tree* t) {
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
    (void)huffman_inc_frequency(t, i);
    // (void) because inability to update frequency is OK here
    return (uint64_t)i;
}

static void squeeze_read_header(bitstream* bs, uint64_t *bytes) {
    uint64_t b = bitstream_read_bits(bs, sizeof(uint64_t) * 8);
    if (bs->error == 0) { *bytes = b; }
}

static uint16_t squeeze_read_length(squeeze_type* s, uint16_t lit) {
    const uint8_t base = (uint8_t)(lit - squeeze_deflate_sym_min);
    if (base >= countof(squeeze_len_base)) {
        s->error = EINVAL;
        return 0;
    } else {
        const uint8_t bits = squeeze_len_xb[base];
        if (bits != 0) {
            uint64_t extra = squeeze_read_bits(s, bits);
            assert(squeeze_len_base[base] + extra <= UINT16_MAX);
            return s->error == 0 ?
                (uint16_t)squeeze_len_base[base] + (uint16_t)extra : 0;
        } else {
            return squeeze_len_base[base];
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
    if (s->error == 0 && base >= countof(squeeze_pos_base)) {
        s->error = EINVAL;
    } else {
        pos = squeeze_pos_base[base];
    }
    if (s->error == 0) {
        uint8_t bits  = squeeze_pos_xb[base];
        if (bits > 0) {
            uint64_t extra = squeeze_read_bits(s, bits);
            if (s->error == 0) { pos += (uint32_t)extra; }
        }
    }
    return pos;
}

static void squeeze_decompress(squeeze_type* s, bitstream* bs,
                               uint8_t* data, uint64_t bytes) {
    s->bs = bs;
    if (!huffman_insert(&s->lit, squeeze_lit_nyt)) { s->error = EINVAL; }
    if (!huffman_insert(&s->pos, squeeze_pos_nyt)) { s->error = EINVAL; }
    squeeze_deflate_init(s);
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
            assert(squeeze_deflate_sym_min <= lit && lit < squeeze_lit_nyt);
            if (squeeze_deflate_sym_min <= lit && lit <= squeeze_lit_nyt) {
                uint32_t len = squeeze_read_length(s, (uint16_t)lit);
                if (s->error != 0) { break; }
                assert(squeeze_deflate_len_min <= len && len <= squeeze_deflate_len_max);
                if (squeeze_deflate_len_min <= len && len <= squeeze_deflate_len_max) {
                    uint32_t pos = squeeze_read_pos(s);
                    if (s->error != 0) { break; }
                    assert(0 < pos && pos <= 0x7FFF);
                    if (0 < pos && pos <= 0x7FFF) {
                        // memcpy() cannot be used on overlapped regions
                        // because it may read more than one byte at a time.
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

#if defined(__cplusplus)
extern "C" { // avoid mangling of global variable "squeeze" by C++ compilers
#endif

squeeze_interface squeeze = {
    .init         = squeeze_init,
    .write_header = squeeze_write_header,
    .compress     = squeeze_compress,
    .read_header  = squeeze_read_header,
    .decompress   = squeeze_decompress,
};

#if defined(__cplusplus)
} // extern "C"
#endif

#endif // squeeze_implementation
