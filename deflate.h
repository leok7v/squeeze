#ifndef deflate_header_included
#define deflate_header_included

#include <stdint.h>

enum {
    deflate_len_count = 285 + 1,
    deflate_pos_count = 32768 + 1
};

static uint8_t deflate_len_index[deflate_len_count];
static uint8_t deflate_pos_index[deflate_pos_count];

static const uint16_t deflate_len_base[29] = {
    3, 4, 5, 6, 7, 8, 9, 10,  // 257-264
    11, 13, 15, 17,           // 265-268
    19, 23, 27, 31,           // 269-272
    35, 43, 51, 59,           // 273-276
    67, 83, 99, 115,          // 277-280
    131, 163, 195, 227, 258   // 281-285
};

static const uint8_t deflate_len_extra_bits[29] = {
    0, 0, 0, 0, 0, 0, 0, 0,   // 257-264
    1, 1, 1, 1,               // 265-268
    2, 2, 2, 2,               // 269-272
    3, 3, 3, 3,               // 273-276
    4, 4, 4, 4,               // 277-280
    5, 5, 5, 5, 0             // 281-285 (len = 258 has no extra bits)
};

static const uint16_t deflate_pos_base[30] = {
    1, 2, 3, 4,               // 0-3
    5, 7,                     // 4-5
    9, 13,                    // 6-7
    17, 25,                   // 8-9
    33, 49,                   // 10-11
    65, 97,                   // 12-13
    129, 193,                 // 14-15
    257, 385,                 // 16-17
    513, 769,                 // 18-19
    1025, 1537,               // 20-21
    2049, 3073,               // 22-23
    4097, 6145,               // 24-25
    8193, 12289,              // 26-27
    16385, 24577              // 28-29
};

static const uint8_t deflate_pos_extra_bits[30] = {
    0, 0, 0, 0,               // 0-3
    1, 1,                     // 4-5
    2, 2,                     // 6-7
    3, 3,                     // 8-9
    4, 4,                     // 10-11
    5, 5,                     // 12-13
    6, 6,                     // 14-15
    7, 7,                     // 16-17
    8, 8,                     // 18-19
    9, 9,                     // 20-21
    10, 10,                   // 22-23
    11, 11,                   // 24-25
    12, 12,                   // 26-27
    13, 13                    // 28-29
};

static void deflate_init(void) {
    uint8_t  j = 0;
    uint16_t n = deflate_len_base[j] + (1u << deflate_len_extra_bits[j]); 
    for (int16_t i = 3; i < deflate_len_count; i++) {
        if (i == n) {
            j++;
            n = deflate_len_base[j] + (1u << deflate_len_extra_bits[j]);
        }
        deflate_len_index[i] = j;
    }
    assert(j == 28);
    j = 0;
    n = deflate_pos_base[j] + (1u << deflate_pos_extra_bits[j]); 
    for (int32_t i = 0; i < deflate_pos_count; i++) {
        if (i == n) {
            j++;
            n = deflate_pos_base[j] + (1u << deflate_pos_extra_bits[j]); 
        }
        deflate_pos_index[i] = j;
    }
    assert(j == 29);
}

/*
    DEFLATE compression (most commonly used compression method in ZIP), 
    the compressed data is represented as a bitstream consisting of a 
    combination of Literals, Length/Distance pairs (for backreferences).

    Overview of the DEFLATE Encoding Bitstream:

    Literals:

    A literal is any single uncompressed byte from the input.
    Literal values are represented by their ASCII code (0-255) 
    and are output directly (via Huffman code) when encountered in the input.
    Length/Distance pairs (Backreferences):

    When a sequence of bytes has been repeated earlier in the input stream, 
    the DEFLATE algorithm uses a backreference to refer to the earlier 
    occurrence.
    This is encoded as a length (for how many bytes to copy) and a distance 
    (where to copy from in the previously decompressed data).

    The bitstream encodes these two values as follows:

    Length:
        The length value is encoded with a base length and extra bits.
        There are 29 possible length codes ranging from 3 to 258.
        Length codes 257-285 are reserved for length encoding. 
        Each length code has a corresponding base length and extra bits 
        to refine the exact length.

    Distance:
        The distance value is encoded similarly, with a base distance 
        and extra bits.
        There are 30 possible distance codes for distances 
        from 1 to 32,768 bytes.
        Distance codes range from 0 to 29, and each code represents a 
        base distance with extra bits to refine the exact distance.
*/

#endif // deflate_header_included
