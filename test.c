#include "rt.h"

#undef assert
#undef swear
#undef countof
#undef min
#undef max
#undef swap

#define assert(b, ...) rt_assert(b, __VA_ARGS__)
#define swear(b, ...)  rt_swear(b, __VA_ARGS__) // release build assert
#define printf(...)    rt_printf(__VA_ARGS__)
#define println(...)   rt_println(__VA_ARGS__)  // same as printf("...\n", ...)
#define countof(a)     rt_countof(a)
#define min(x, y)      rt_min(x, y)
#define max(x, y)      rt_max(x, y)
#define swap(a, b)     rt_swap(a, b)

#include "bitstream.h"
#include "map.h"
#include "squeeze.h"
#include "file.h"

// enum { bits_win = 12, bits_map = 19, bits_len = 4 };

//enum { bits_win = 12, bits_map = 16, bits_len = 8 };

enum { bits_win = 15, bits_map = 19, bits_len = 8 };

static squeeze_type* squeeze_new(bitstream_type* bs, uint8_t win_bits,
                                 uint8_t map_bits, uint8_t len_bits) {
    const uint64_t bytes = squeeze_sizeof(win_bits, map_bits, len_bits);
    squeeze_type* s = (squeeze_type*)calloc(1, (size_t)bytes);
    if (s != null) {
        squeeze.init_with(s, s, bytes, win_bits, map_bits, len_bits);
        s->bs = bs;
    }
    return s;
}

static void squeeze_delete(squeeze_type* s) {
    free(s);
}

static errno_t compress(const char* from, const char* to,
                        const uint8_t* data, uint64_t bytes) {
    FILE* out = null; // compressed file
    errno_t r = fopen_s(&out, to, "wb") != 0;
    if (r != 0 || out == null) {
        printf("Failed to create \"%s\": %s\n", to, strerror(r));
        return r;
    }
    squeeze_type* s = null;
    bitstream_type bs = { .file = out };
    squeeze.write_header(&bs, bytes, bits_win, bits_map, bits_len);
    if (bs.error != 0) {
        r = bs.error;
        printf("Failed to create \"%s\": %s\n", to, strerror(r));
    } else {
        s = squeeze_new(&bs, bits_win, bits_map, bits_len);
        if (s != null) {
            squeeze.compress(s, data, bytes, 1u << bits_win);
            assert(s->error == 0);
        } else {
            r = ENOMEM;
            printf("squeeze_new() failed.\n");
            assert(false);
        }
    }
    errno_t rc = fclose(out) == 0 ? 0 : errno; // error writing buffered output
    if (rc != 0) {
        printf("Failed to flush on file close: %s\n", strerror(rc));
        if (r == 0) { r = rc; }
    }
    if (r == 0) {
        r = s->error;
        if (r != 0) {
            printf("Failed to compress: %s\n", strerror(r));
        } else {
            char* fn = from == null ? null : strrchr(from, '\\'); // basename
            if (fn == null) { fn = from == null ? null : strrchr(from, '/'); }
            if (fn != null) { fn++; } else { fn = (char*)from; }
            const uint64_t written = s->bs->bytes;
            double percent = written * 100.0 / bytes;
            if (from != null) {
                printf("%7lld -> %7lld %5.1f%% of \"%s\"\n", bytes, written,
                                                             percent, fn);
            } else {
                printf("%7lld -> %7lld %5.1f%%\n", bytes, written, percent);
            }
        }
    }
    if (s != null) {
        squeeze_delete(s); s = null;
    }
    return r;
}

static errno_t verify(const char* fn, const uint8_t* input, size_t size) {
    // decompress and compare
    FILE* in = null; // compressed file
    errno_t r = fopen_s(&in, fn, "rb");
    if (r != 0 || in == null) {
        printf("Failed to open \"%s\"\n", fn);
    }
    bitstream_type bs = { .file = in };
    uint64_t bytes = 0;
    uint8_t win_bits = 0;
    uint8_t map_bits = 0;
    uint8_t len_bits = 0;
    if (r == 0) {
        squeeze.read_header(&bs, &bytes, &win_bits, &map_bits, &len_bits);
        if (bs.error != 0) {
            printf("Failed to read header from \"%s\"\n", fn);
            r = bs.error;
        }
    }
    if (r == 0) {
        squeeze_type* s = squeeze_new(&bs, win_bits, map_bits, len_bits);
        if (s == null) {
            r = ENOMEM;
            printf("squeeze_new() failed.\n");
            assert(false);
        } else {
            assert(s->error == 0 && bytes == size && win_bits == win_bits);
            uint8_t* data = (uint8_t*)calloc(1, (size_t)bytes);
            if (data == null) {
                printf("Failed to allocate memory for decompressed data\n");
                fclose(in);
                return ENOMEM;
            }
            squeeze.decompress(s, data, bytes);
            fclose(in);
            assert(s->error == 0);
            if (s->error == 0) {
                const bool same = size == bytes && memcmp(input, data, bytes) == 0;
                assert(same);
                if (!same) {
                    printf("compress() and decompress() are not the same\n");
                    // ENODATA is not original posix error but is OpenGroup error
                    r = ENODATA; // or EIO
                } else if (bytes < 128) {
                    printf("decompressed: %.*s\n", (unsigned int)bytes, data);
                }
            }
            free(data);
            if (r != 0) {
                printf("Failed to decompress\n");
            }
            squeeze_delete(s); s = null;
        }
    }
    return r;
}

const char* compressed = "~compressed~.bin";

static errno_t test(const char* fn, const uint8_t* data, size_t bytes) {
    errno_t r = compress(fn, compressed, data, bytes);
    if (r == 0) {
        r = verify(compressed, data, bytes);
    }
    (void)remove(compressed);
    return r;
}

static errno_t test_compression(const char* fn) {
    uint8_t* data = null;
    size_t bytes = 0;
    errno_t r = file.read_fully(fn, &data, &bytes);
    if (r != 0) { return r; }
    return test(fn, data, bytes);
}

static errno_t locate_test_folder(void) {
    // on Unix systems with "make" executable usually resided
    // and is run from root of repository... On Windows with
    // MSVC it is buried inside bin/... folder depths
    // on X Code in MacOS it can be completely out of tree.
    // So we need to find the test files.
    for (;;) {
        if (file.exist("test/bible.txt")) { return 0; }
        if (file.chdir("..") != 0) { return errno; }
    }
}

int main(int argc, const char* argv[]) {
    (void)argc; (void)argv; // unused
    errno_t r = locate_test_folder();
    if (r == 0) {
        uint8_t data[4 * 1024] = {0};
        r = test(null, data, sizeof(data));
        // lz77 deals with run length encoding in amazing overlapped way
        for (int32_t i = 0; i < sizeof(data); i += 4) {
            memcpy(data + i, "\x01\x02\x03\x04", 4);
        }
        r = test(null, data, sizeof(data));
    }
    if (r == 0) {
        const char* data = "Hello World Hello.World Hello World";
        size_t bytes = strlen((const char*)data);
        r = test(null, (const uint8_t*)data, bytes);
    }
    if (r == 0 && file.exist(__FILE__)) { // test.c source code:
        r = test_compression(__FILE__);
    }
    // argv[0] executable filepath (Windows) or possibly name (Unix)
    if (r == 0 && file.exist(argv[0])) {
        r = test_compression(argv[0]);
    }
    static const char* test_files[] = {
        "test/bible.txt",     // bits len:3.01 pos:10.73 #words:91320 #lens:112
        "test/hhgttg.txt",    // bits len:2.33 pos:10.78 #words:12034 #lens:40
        "test/confucius.txt",
        "test/laozi.txt",
        "test/sqlite3.c",
        "test/arm64.elf",
        "test/x64.elf",
        "test/mandrill.bmp",
        "test/mandrill.png",
    };
    for (int i = 0; i < countof(test_files) && r == 0; i++) {
        if (file.exist(test_files[i])) {
            r = test_compression(test_files[i]);
        }
    }
    return r;
}

#define map_implementation
#include "map.h"

#define bitstream_implementation
#include "bitstream.h"

#define huffman_implementation
#include "huffman.h"

#define file_implementation
#include "file.h"

#define squeeze_implementation
#include "squeeze.h"

#if 0

WITHOUT HUFFMAN:

     35 ->      32  91.4%
   4096 ->      16   0.4%
   4096 ->      24   0.6%
  10502 ->    3928  37.4% of "test.c"
 229376 ->  129248  56.3% of "sqz.exe"
4436173 -> 2334152  52.6% of "bible.txt"
 279056 ->  179992  64.5% of "hhgttg.txt"
  68762 ->   39232  57.1% of "confucius.txt"
  20943 ->   12728  60.8% of "laozi.txt"
8182289 -> 3519736  43.0% of "sqlite3.c"
 847400 ->  522224  61.6% of "arm64.elf"
 926536 ->  562440  60.7% of "x64.elf"
 786570 ->  830336 105.6% of "mandrill.bmp"
 627896 ->  667472 106.3% of "mandrill.png"

WITH HUFFMAN + DICTIONARY:

     35 ->      40 114.3%
   4096 ->      24   0.6%
   4096 ->      24   0.6%
  10716 ->    3616  33.7% of "test.c"
 231424 ->  116680  50.4% of "sqz.exe"
4436173 -> 1801336  40.6% of "bible.txt"
 279056 ->  140952  50.5% of "hhgttg.txt"
  68762 ->   32424  47.2% of "confucius.txt"
  20943 ->   10664  50.9% of "laozi.txt"
8182289 -> 2707072  33.1% of "sqlite3.c"
 847400 ->  456024  53.8% of "arm64.elf"
 926536 ->  514016  55.5% of "x64.elf"
 786570 ->  910648 115.8% of "mandrill.bmp"
 627896 ->  747184 119.0% of "mandrill.png"

// bible.txt.zip size (zip: 1,398,871 bytes)
// enum { bits_win = 12, bits_map = 19, bits_len = 4 }; // 1801336
// enum { bits_win = 12, bits_map = 14, bits_len = 5 }; // 1736248
// enum { bits_win = 13, bits_map = 15, bits_len = 6 }; // 1683192
// enum { bits_win = 14, bits_map = 16, bits_len = 7 }; // 1623176
// enum { bits_win = 15, bits_map = 17, bits_len = 8 }; // 1559648
// enum { bits_win = 16, bits_map = 18, bits_len = 8 }; // 1496552


#endif
