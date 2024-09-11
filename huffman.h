#ifndef huffman_header_included
#define huffman_header_included

#include <stdint.h>
#include <errno.h>

// Adaptive Huffman Coding
// https://en.wikipedia.org/wiki/Adaptive_Huffman_coding

typedef struct huffman_node_struct {
    uint64_t freq;
    uint64_t path;
    int32_t  bits; // 0 for root
    int32_t  pix;  // parent
    int32_t  lix;  // left
    int32_t  rix;  // right
} huffman_node_type;

typedef struct huffman_tree_struct {
    huffman_node_type* node;
    int32_t n;
    int32_t depth; // max tree depth seen
    int32_t complete; // tree is too deep or freq too high - no more updates
    int32_t padding;
    // stats:
    struct {
        size_t updates;
        size_t swaps;
        size_t moves;
    } stats;
} huffman_tree_type;

typedef struct {
    void (*init)(huffman_tree_type* t, huffman_node_type nodes[],
                 const size_t m);
    void (*inc_frequency)(huffman_tree_type* t, int32_t symbol);
    double (*entropy)(const huffman_tree_type* t); // Shannon entropy (bps)
    uint8_t (*log2_of_pow2)(uint64_t pow2);
} huffman_interface;

extern huffman_interface huffman;

#endif // huffman_header_included

#if defined(huffman_implementation) && !defined(huffman_implemented)

#define huffman_implemented

#ifndef assert
#include <assert.h>
#endif

#ifndef null
#define null ((void*)0)
#endif

static void huffman_update_paths(huffman_tree_type* t, int32_t i) {
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
        assert(rix != -1);
        t->node[lix].bits = bits + 1;
        t->node[lix].path = path;
        t->node[rix].bits = bits + 1;
        t->node[rix].path = path | (1ULL << bits);
        huffman_update_paths(t, lix);
        huffman_update_paths(t, rix);
    } else {
        if (bits > t->depth) { t->depth = bits; }
    }
}

static inline int32_t huffman_swap_siblings(huffman_tree_type* t,
                                            const int32_t ix) {
    const int32_t m = t->n * 2 - 1;
    assert(0 <= ix && ix < m);
    if (ix < m - 1) { // not root
        const int32_t pix = t->node[ix].pix; // parent (cannot be a leaf)
        const int32_t lix = t->node[pix].lix; assert(0 <= lix && lix < m - 1);
        const int32_t rix = t->node[pix].rix; assert(0 <= rix && rix < m - 1);
        if (t->node[lix].freq > t->node[rix].freq) { // swap
            t->stats.swaps++;
            t->node[pix].lix = rix;
            t->node[pix].rix = lix;
            huffman_update_paths(t, pix); // because swap changed all path below
            return ix == lix ? rix : lix;
        }
    }
    return ix;
}

static inline void huffman_frequency_changed(huffman_tree_type* t, int32_t ix);

static inline void huffman_update_freq(huffman_tree_type* t, int32_t ix) {
    const int32_t lix = t->node[ix].lix; assert(lix != -1);
    const int32_t rix = t->node[ix].rix; assert(rix != -1);
    t->node[ix].freq = t->node[lix].freq + t->node[rix].freq;
}

static inline void huffman_move_up(huffman_tree_type* t, int32_t ix) {
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

static void huffman_frequency_changed(huffman_tree_type* t, int32_t i) {
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

static inline void huffman_inc_frequency(huffman_tree_type* t, int32_t ix) {
    assert(0 <= ix && ix < t->n); // terminal
    // If input sequence frequencies are severely skewed (e.g. Lucas numbers
    // similar to Fibonacci numbers) and input sequence is long enough.
    // The depth of the tree will grow past 64 bits.
    // The first Lucas number that exceeds 2^64 is
    // L(81) = 18,446,744,073,709,551,616 not actually realistic but
    // better be safe than sorry:
    if (!t->complete) {
        if (t->depth < 63 && t->node[ix].freq < UINT64_MAX - 1) {
            t->node[ix].freq++;
            huffman_frequency_changed(t, ix);
        } else {
            // ignore future frequency updates
            t->complete = 1;
        }
    }
}

static uint8_t huffman_log2_of_pow2(uint64_t pow2) {
    assert(pow2 > 0 && (pow2 & (pow2 - 1)) == 0);
    if (pow2 > 0 && (pow2 & (pow2 - 1)) == 0) {
        uint8_t bit = 0;
        while (pow2 >>= 1) { bit++; }
        return bit;
    } else {
        return 0xFF; // error
    }
}

static double huffman_entropy(const huffman_tree_type* t) { // Shannon entropy
    double total = 0;
    double aha_entropy = 0.0;
    for (int32_t i = 0; i < t->n; i++) { total += (double)t->node[i].freq; }
    for (int32_t i = 0; i < t->n; i++) {
        if (t->node[i].freq > 0) {
            double p_i = (double)t->node[i].freq / total;
            aha_entropy += p_i * log2(p_i);
        }
    }
    return -aha_entropy;
}

static void huffman_init(huffman_tree_type* t, huffman_node_type nodes[],
                         const size_t count) {
    assert(7 <= count && count < INT32_MAX); // must pow(2, bits_per_symbol) * 2 - 1
    const int32_t m = (int32_t)count;
    const int32_t n = (int32_t)(m + 1) / 2;
    assert(n > 4 && (n & (n - 1)) == 0); // must be power of 2
    const int32_t bits_per_symbol = huffman_log2_of_pow2(n);
    assert(2 <= bits_per_symbol && bits_per_symbol <= 20);
    memset(&t->stats, 0x00, sizeof(t->stats));
    t->node = nodes;
    t->n = n;
    t->depth = bits_per_symbol;
    t->complete = 0;
    for (int32_t i = 0; i < n; i++) {
        t->node[i] = (huffman_node_type){
            .freq = 1, .lix = -1, .rix = -1, .pix = n + i / 2,
            .bits = bits_per_symbol
        };
    }
    int32_t ix = n;
    int32_t lix = 0;
    int32_t rix = 1;
    int32_t n2  = n / 2;
    int32_t bits = bits_per_symbol - 1;
    while (n2 > 0) {
        int32_t pix = ix + n2;
        for (int32_t i = 0; i < n2; i++) {
            uint64_t f = t->node[lix].freq + t->node[rix].freq;
            assert(ix < m);
            t->node[ix] = (huffman_node_type){
                .freq = f, .lix = lix, .rix = rix, .pix = pix, .bits = bits };
            lix += 2;
            rix += 2;
            if (i % 2 == 1) { pix++; }
            ix++;
        }
        n2 = n2 / 2;
        bits--;
    }
    // change root parent to be -1
    const int32_t root = m - 1;
    assert(t->node[root].bits == 0);
    assert(t->node[root].pix == m);
    t->node[root].pix = -1;
    t->node[root].path = 0;
    huffman_update_paths(t, m - 1);
}

huffman_interface huffman = {
    .init          = huffman_init,
    .inc_frequency = huffman_inc_frequency,
    .log2_of_pow2  = huffman_log2_of_pow2,
    .entropy       = huffman_entropy
};

#endif