#ifndef file_header_included
#define file_header_included

#include <errno.h>
#include <stdint.h>

typedef struct {
    errno_t (*chdir)(const char* name);
    bool    (*exist)(const char* filename);
    errno_t (*size)(FILE* f, size_t* size);
    errno_t (*read_fully)(const char* fn, const uint8_t* *data, size_t *bytes);
} file_interface;

extern file_interface file;

#endif // file_header_included

#if defined(file_implementation) && !defined(file_implemented)

#define file_implemented

#include <stdio.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h> // chdir
#else
#include <unistd.h> // chdir
#endif

#ifndef null
#define null ((void*)0)
#endif

#ifndef assert
#include <assert.h>
#endif

static bool file_exist(const char* filename) {
    struct stat st = {0};
    return stat(filename, &st) == 0;
}

static errno_t file_size(FILE* f, size_t* size) {
    // on error returns (fpos_t)-1 and sets errno
    fpos_t pos = 0;
    if (fgetpos(f, &pos) != 0) { return errno; }
    if (fseek(f, 0, SEEK_END) != 0) { return errno; }
    fpos_t eof = 0;
    if (fgetpos(f, &eof) != 0) { return errno; }
    if (fseek(f, 0, SEEK_SET) != 0) { return errno; }
    if ((uint64_t)eof > SIZE_MAX) { return E2BIG; }
    *size = (size_t)eof;
    return 0;
}

static errno_t read_whole_file(FILE* f, const uint8_t* *data, size_t *bytes) {
    size_t size = 0;
    errno_t r = file_size(f, &size);
    if (r != 0) { return r; }
    if (size > SIZE_MAX) { return E2BIG; }
    uint8_t* p = (uint8_t*)malloc(size); // does set errno on failure
    if (p == null) { return errno; }
    if (fread(p, 1, size, f) != (size_t)size) { free(p); return errno; }
    *data = p;
    *bytes = (size_t)size;
    return 0;
}

static errno_t read_fully(const char* fn, const uint8_t* *data, size_t *bytes) {
    FILE* f = null;
    errno_t r = fopen_s(&f, fn, "rb");
    if (r != 0) {
        printf("Failed to open file \"%s\": %s\n", fn, strerror(r));
        return r;
    }
    r = read_whole_file(f, data, bytes); // to the heap
    if (r != 0) {
        printf("Failed to read file \"%s\": %s\n", fn, strerror(r));
        fclose(f);
        return r;
    }
    return fclose(f) == 0 ? 0 : errno;
}

static errno_t file_chdir(const char* name) {
    if (chdir(name) != 0) { return errno; }
    return 0;
}

file_interface file = {
    .chdir      = file_chdir,
    .exist      = file_exist,
    .size       = file_size,
    .read_fully = read_fully
};

#endif // file_implementation
