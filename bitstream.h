#ifndef bitstream_header_included
#define bitstream_header_included

#include <errno.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

typedef struct bitstream_struct {
    FILE*    file; // file and (data,capacity) is exclusive
    uint8_t* data;
    uint64_t capacity; // data[capacity]
    uint64_t bytes; // number of bytes written
    uint64_t read;  // number of bytes read
    uint64_t b64;   // bit shifting buffer
    int32_t  bits;  // bit count inside b64
    errno_t  error; // sticky error
} bitstream_type;

static inline void     bitstream_create(bitstream_type* bs, void* data, size_t capacity);
static inline void     bitstream_write_bit(bitstream_type* bs, int32_t bit);
static inline void     bitstream_write_bits(bitstream_type* bs, uint64_t data, int32_t bits);
static inline bool     bitstream_read_bit(bitstream_type* bs);
static inline uint64_t bitstream_read_bits(bitstream_type* bs, int32_t bits);
static inline void     bitstream_flush(bitstream_type* bs); // write trailing zeros
static inline void     bitstream_dispose(bitstream_type* bs);
#endif // bitstream_header_included

#if defined(bitstream_implementation) && !defined(bitstream_implemented)

#define bitstream_implemented

static inline void bitstream_write_bit(bitstream_type* bs, int32_t bit) {
    if (bs->error == 0) {
        bs->b64 <<= 1;
        bs->b64 |= (bit & 1);
        bs->bits++;
        if (bs->bits == 64) {
            if (bs->data != null && bs->capacity > 0) {
                assert(bs->file == null);
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
                size_t written = fwrite(&bs->b64, 1, 8, bs->file);
                bs->error = written == 8 ? 0 : errno;
                if (bs->error == 0) { bs->bytes += 8; }
            }
            bs->bits = 0;
            bs->b64 = 0;
        }
    }
}

static inline void bitstream_write_bits(bitstream_type* bs, 
                                        uint64_t data, int32_t bits) {
    assert(0 < bits && bits <= 64);
    while (bits > 0 && bs->error == 0) {
        bitstream_write_bit(bs, data & 1);
        bits--;
        data >>= 1;
    }
}

static inline bool bitstream_read_bit(bitstream_type* bs) {
    bool bit = false;
    if (bs->error == 0) {
        if (bs->bits == 0) {
            bs->b64 = 0;
            if (bs->data != null && bs->bytes > 0) {
                assert(bs->file == null);
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
                size_t read = fread(&bs->b64, 1, 8, bs->file);
                bs->error = read == 8 ? 0 : errno;
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

static inline uint64_t bitstream_read_bits(bitstream_type* bs, int32_t bits) {
    uint64_t data = 0;
    assert(0 < bits && bits <= 64);
    for (int32_t b = 0; b < bits && bs->error == 0; b++) {
        bool bit = bitstream_read_bit(bs);
        if (bit) { data |= ((uint64_t)bit) << b; }
    }
    return data;
}

static inline void bitstream_create(bitstream_type* bs, void* data, size_t capacity) {
    assert(bs->data != null);
    memset(bs, 0x00, sizeof(*bs));
    bs->data = (uint8_t*)data;
    bs->capacity  = capacity;
}

static inline void bitstream_flush(bitstream_type* bs) {
    while (bs->bits > 0 && bs->error == 0) { bitstream_write_bit(bs, 0); }
}

static inline void bitstream_dispose(bitstream_type* bs) {
    memset(bs, 0x00, sizeof(*bs));
}

#endif // bitstream_implementation