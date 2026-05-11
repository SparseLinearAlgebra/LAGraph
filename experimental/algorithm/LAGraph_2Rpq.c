//
// Different RPQ semantics
//

#define LG_FREE_WORK                            \
{                                               \
}
#define LG_FREE_ALL                         \
{                                           \
        LG_FREE_WORK ;                          \
}

#include "LG_internal.h"
#include "LAGraphX.h"

#include <assert.h>
#include <limits.h>

#include <stdatomic.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PATH_LIMIT 100000

#define RPQ_MAX_PATH_LENGTH 60000

#define RPQ_STRINGIFY_HELPER(x) #x
#define RPQ_STRINGIFY(x) RPQ_STRINGIFY_HELPER(x)

// This define and three functions below need for YAGO dataset.
// Because we have OOM on it
// 
#define PATHS_PER_POINT_LIMIT 3000000

static atomic_int path_limit_exceeded ;

static void path_limit_reset (void)
{
    atomic_store (&path_limit_exceeded, 0) ;
}

static void path_limit_mark_exceeded (void)
{
    atomic_store (&path_limit_exceeded, 1) ;
}

static bool path_limit_failed (void)
{
    return atomic_load (&path_limit_exceeded) != 0 ;
}
//

// Temporary storage for data referenced by GraphBLAS UDT values.
// GraphBLAS copies MultiplePaths by value, so pointers inside it must point to
// memory that remains alive while frontier matrices are alive.
#define TEMP_ARENA_BYTES ((size_t) 32 * 1024 * 1024 * 1024)

static unsigned char *temp_arena_data = NULL ;
static size_t temp_arena_capacity = 0 ;
static atomic_size_t temp_arena_offset ;
static atomic_int temp_arena_oom ;

static size_t align_up_size (size_t value, size_t align)
{
    size_t rem = value % align ;

    if (rem == 0)
    {
        return value ;
    }

    return value + (align - rem) ;
}

static int temp_arena_init (char *msg)
{
    int info ;

    temp_arena_data = NULL ;
    temp_arena_capacity = 0 ;
    atomic_store (&temp_arena_offset, 0) ;
    atomic_store (&temp_arena_oom, 0) ;

    info = LAGraph_Malloc ((void **) &temp_arena_data,
        TEMP_ARENA_BYTES, sizeof (unsigned char), msg) ;
    if (info != GrB_SUCCESS)
    {
        atomic_store (&temp_arena_oom, 1) ;
        return info ;
    }

    temp_arena_capacity = TEMP_ARENA_BYTES ;
    return GrB_SUCCESS ;
}

static void temp_arena_destroy (void)
{
    LAGraph_Free ((void **) &temp_arena_data, NULL) ;
    temp_arena_capacity = 0 ;
    atomic_store (&temp_arena_offset, 0) ;
    atomic_store (&temp_arena_oom, 0) ;
}

static void *temp_alloc (size_t size)
{
    size_t align ;
    size_t aligned_size ;
    size_t start ;

    if (size == 0)
    {
        return NULL ;
    }

    align = _Alignof (max_align_t) ;
    aligned_size = align_up_size (size, align) ;
    start = atomic_fetch_add (&temp_arena_offset, aligned_size) ;

    if (start > temp_arena_capacity ||
        aligned_size > temp_arena_capacity - start)
    {
        atomic_store (&temp_arena_oom, 1) ;
        return NULL ;
    }

    return (void *) (temp_arena_data + start) ;
}

static void *temp_calloc_bytes (size_t size)
{
    void *ptr ;

    ptr = temp_alloc (size) ;
    if (ptr != NULL)
    {
        memset (ptr, 0, size) ;
    }

    return ptr ;
}

static bool temp_alloc_failed (void)
{
    return atomic_load (&temp_arena_oom) != 0 ;
}

static size_t temp_arena_used (void)
{
    return atomic_load (&temp_arena_offset) ;
}
//



// JIT kernels are compiled into a separate shared object. They cannot call
// static functions from this translation unit, so expose tiny wrappers for
// the stateful parts: arena allocation and path-limit reporting.
LAGRAPHX_PUBLIC
void *LAGraph_Rpq_jit_temp_calloc_bytes (size_t size)
{
    return temp_calloc_bytes (size) ;
}

LAGRAPHX_PUBLIC
void LAGraph_Rpq_jit_path_limit_mark_exceeded (void)
{
    path_limit_mark_exceeded () ;
}
//

typedef struct MultiplePathsExtra {
    size_t count;
    Path paths[0];
} MultiplePathsExtra;

typedef struct {
    Path paths[QUICK_PATH_COUNT];
    size_t path_count;
    MultiplePathsExtra *extra;
} MultiplePaths ;

MultiplePaths multiple_paths_identity ;

GrB_Type multiple_paths ;
GrB_BinaryOp combine_multiple_paths_op ;
GrB_Monoid combine_multiple_paths ;
GrB_BinaryOp first_multiple_paths ;
GrB_BinaryOp second_multiple_paths ;
GrB_Semiring first_combine_multiple_paths ;
GrB_Semiring second_combine_multiple_paths ;
GrB_IndexUnaryOp extend_multiple_paths ;
GrB_IndexUnaryOp extend_multiple_simple ;
GrB_IndexUnaryOp extend_multiple_trails ;




// JIT definitions
#define MULTIPLE_PATHS_TYPE_DEFN                                      \
"#include <stdint.h>\n"                                                \
"#include <stddef.h>\n"                                                \
"#include <stdbool.h>\n"                                               \
"#include <string.h>\n"                                                \
"#define QUICK_PATH_LENGTH " RPQ_STRINGIFY(QUICK_PATH_LENGTH) "\n"     \
"#define QUICK_PATH_COUNT " RPQ_STRINGIFY(QUICK_PATH_COUNT) "\n"       \
"#define PATHS_PER_POINT_LIMIT " RPQ_STRINGIFY(PATHS_PER_POINT_LIMIT) "\n"\
"#define RPQ_MAX_PATH_LENGTH " RPQ_STRINGIFY(RPQ_MAX_PATH_LENGTH) "\n" \
"typedef uint64_t Vertex;\n"                                           \
"typedef struct PathExtra {\n"                                         \
"    size_t len;\n"                                                    \
"    Vertex vertices[];\n"                                             \
"} PathExtra;\n"                                                       \
"typedef struct Path {\n"                                              \
"    Vertex vertices[QUICK_PATH_LENGTH];\n"                            \
"    size_t vertex_count;\n"                                           \
"    PathExtra *extra;\n"                                              \
"} Path;\n"                                                           \
"typedef struct MultiplePathsExtra {\n"                                \
"    size_t count;\n"                                                  \
"    Path paths[];\n"                                                  \
"} MultiplePathsExtra;\n"                                              \
"typedef struct MultiplePaths {\n"                                     \
"    Path paths[QUICK_PATH_COUNT];\n"                                  \
"    size_t path_count;\n"                                             \
"    MultiplePathsExtra *extra;\n"                                     \
"} MultiplePaths;\n"                                                   \
"extern void *LAGraph_Rpq_jit_temp_calloc_bytes(size_t size);\n"        \
"extern void LAGraph_Rpq_jit_path_limit_mark_exceeded(void);\n"         \
"static size_t path_extra_len_jit(const Path *path)\n"                 \
"{\n"                                                                  \
"    return path->vertex_count > QUICK_PATH_LENGTH ?\n"                \
"        path->vertex_count - QUICK_PATH_LENGTH : 0;\n"                \
"}\n"                                                                  \
"static Vertex path_get_vertex_jit(const Path *path, size_t i)\n"       \
"{\n"                                                                  \
"    if (i < QUICK_PATH_LENGTH) return path->vertices[i];\n"            \
"    return path->extra->vertices[i - QUICK_PATH_LENGTH];\n"           \
"}\n"                                                                  \
"static bool path_is_empty_jit(const Path *path)\n"                    \
"{\n"                                                                  \
"    return path->vertex_count == 0;\n"                                 \
"}\n"                                                                  \
"static Vertex path_start_vertex_jit(const Path *path)\n"              \
"{\n"                                                                  \
"    return path_get_vertex_jit(path, 0);\n"                            \
"}\n"                                                                  \
"static Vertex path_last_vertex_jit(const Path *path)\n"               \
"{\n"                                                                  \
"    return path_get_vertex_jit(path, path->vertex_count - 1);\n"       \
"}\n"                                                                  \
"static bool path_is_closed_cycle_jit(const Path *path)\n"             \
"{\n"                                                                  \
"    return path->vertex_count > 1 &&\n"                                \
"        path_start_vertex_jit(path) == path_last_vertex_jit(path);\n"  \
"}\n"                                                                  \
"static bool all_paths_can_extend_jit(const Path *path)\n"             \
"{\n"                                                                  \
"    return RPQ_MAX_PATH_LENGTH == 0 ||\n"                              \
"        path->vertex_count < RPQ_MAX_PATH_LENGTH;\n"                   \
"}\n"                                                                  \
"static PathExtra *path_extra_temp_alloc_copy_plus_one_jit(\n"          \
"    const Path *src, Vertex vertex)\n"                                 \
"{\n"                                                                  \
"    size_t old_extra_len = path_extra_len_jit(src);\n"                 \
"    size_t new_extra_len = old_extra_len + 1;\n"                       \
"    size_t nbytes = sizeof(PathExtra) +\n"                             \
"        new_extra_len * sizeof(Vertex);\n"                             \
"    PathExtra *extra = (PathExtra *)\n"                                \
"        LAGraph_Rpq_jit_temp_calloc_bytes(nbytes);\n"                  \
"    if (extra == NULL) return NULL;\n"                                  \
"    extra->len = new_extra_len;\n"                                     \
"    if (old_extra_len > 0)\n"                                          \
"    {\n"                                                              \
"        memcpy(extra->vertices, src->extra->vertices,\n"               \
"            old_extra_len * sizeof(Vertex));\n"                        \
"    }\n"                                                              \
"    extra->vertices[old_extra_len] = vertex;\n"                        \
"    return extra;\n"                                                   \
"}\n"                                                                  \
"static void path_extend_jit(Path *path, Vertex vertex)\n"              \
"{\n"                                                                  \
"    if (path->vertex_count == 0) return;\n"                            \
"    if (path->vertex_count < QUICK_PATH_LENGTH)\n"                     \
"    {\n"                                                              \
"        path->vertices[path->vertex_count] = vertex;\n"                \
"        path->vertex_count++;\n"                                       \
"        return;\n"                                                     \
"    }\n"                                                              \
"    PathExtra *new_extra =\n"                                          \
"        path_extra_temp_alloc_copy_plus_one_jit(path, vertex);\n"      \
"    if (new_extra == NULL)\n"                                          \
"    {\n"                                                              \
"        path->vertex_count = 0;\n"                                     \
"        path->extra = NULL;\n"                                         \
"        return;\n"                                                     \
"    }\n"                                                              \
"    path->extra = new_extra;\n"                                        \
"    path->vertex_count++;\n"                                           \
"}\n"                                                                  \
"static size_t multiple_paths_capacity_jit(const MultiplePaths *x)\n"   \
"{\n"                                                                  \
"    size_t cap = QUICK_PATH_COUNT;\n"                                  \
"    if (x->extra != NULL) cap += x->extra->count;\n"                  \
"    return cap;\n"                                                     \
"}\n"                                                                  \
"static const Path *multiple_paths_nth_const_jit(\n"                   \
"    const MultiplePaths *x, size_t j)\n"                               \
"{\n"                                                                  \
"    if (j < QUICK_PATH_COUNT) return &x->paths[j];\n"                  \
"    return &x->extra->paths[j - QUICK_PATH_COUNT];\n"                  \
"}\n"                                                                  \
"static Path *multiple_paths_nth_mut_jit(MultiplePaths *x, size_t j)\n" \
"{\n"                                                                  \
"    if (j < QUICK_PATH_COUNT) return &x->paths[j];\n"                  \
"    return &x->extra->paths[j - QUICK_PATH_COUNT];\n"                  \
"}\n"                                                                  \
"static MultiplePathsExtra *multiple_paths_extra_temp_alloc_jit(\n"     \
"    size_t extra_count)\n"                                             \
"{\n"                                                                  \
"    if (extra_count == 0) return NULL;\n"                              \
"    size_t nbytes = sizeof(MultiplePathsExtra) +\n"                    \
"        extra_count * sizeof(Path);\n"                                  \
"    MultiplePathsExtra *extra = (MultiplePathsExtra *)\n"              \
"        LAGraph_Rpq_jit_temp_calloc_bytes(nbytes);\n"                  \
"    if (extra == NULL) return NULL;\n"                                  \
"    extra->count = extra_count;\n"                                     \
"    return extra;\n"                                                   \
"}\n"                                                                  \
"static bool multiple_paths_prepare_jit(MultiplePaths *x,\n"            \
"    size_t capacity)\n"                                                \
"{\n"                                                                  \
"    memset(x, 0, sizeof(*x));\n"                                       \
"    if (capacity > PATHS_PER_POINT_LIMIT)\n"                           \
"    {\n"                                                              \
"        LAGraph_Rpq_jit_path_limit_mark_exceeded();\n"                 \
"        return false;\n"                                               \
"    }\n"                                                              \
"    if (capacity > QUICK_PATH_COUNT)\n"                                \
"    {\n"                                                              \
"        x->extra = multiple_paths_extra_temp_alloc_jit(\n"             \
"            capacity - QUICK_PATH_COUNT);\n"                           \
"        if (x->extra == NULL)\n"                                       \
"        {\n"                                                          \
"            x->path_count = 0;\n"                                      \
"            return false;\n"                                           \
"        }\n"                                                          \
"    }\n"                                                              \
"    return true;\n"                                                    \
"}\n"                                                                  \
"static bool multiple_paths_append_jit(MultiplePaths *x,\n"             \
"    const Path *path)\n"                                               \
"{\n"                                                                  \
"    if (x->path_count >= PATHS_PER_POINT_LIMIT)\n"                     \
"    {\n"                                                              \
"        LAGraph_Rpq_jit_path_limit_mark_exceeded();\n"                 \
"        return false;\n"                                               \
"    }\n"                                                              \
"    if (x->path_count >= multiple_paths_capacity_jit(x))\n"            \
"    {\n"                                                              \
"        return false;\n"                                               \
"    }\n"                                                              \
"    *multiple_paths_nth_mut_jit(x, x->path_count) = *path;\n"          \
"    x->path_count++;\n"                                                \
"    return true;\n"                                                    \
"}\n"                                                                  \
"static void multiple_paths_append_unchecked_jit(MultiplePaths *x,\n"   \
"    const Path *path)\n"                                               \
"{\n"                                                                  \
"    *multiple_paths_nth_mut_jit(x, x->path_count) = *path;\n"          \
"    x->path_count++;\n"                                                \
"}\n"                                                                  \
"static bool path_extending_will_add_repeated_vertex_jit(\n"            \
"    const Path *path, Vertex vertex)\n"                                \
"{\n"                                                                  \
"    if (path->vertex_count == 0) return false;\n"                      \
"    if (path_is_closed_cycle_jit(path)) return true;\n"                \
"    for (size_t i = 0 ; i < path->vertex_count ; i++)\n"                 \
"    {\n"                                                              \
"        if (path_get_vertex_jit(path, i) == vertex) return true;\n"    \
"    }\n"                                                              \
"    return false;\n"                                                   \
"}\n"                                                                  \
"static bool path_extending_will_add_repeated_edge_jit(\n"              \
"    const Path *path, Vertex vertex_2)\n"                              \
"{\n"                                                                  \
"    if (path->vertex_count == 0) return false;\n"                      \
"    Vertex vertex_1 = path_last_vertex_jit(path);\n"                   \
"    for (size_t i = 0; i + 1 < path->vertex_count; i++)\n"             \
"    {\n"                                                              \
"        if (path_get_vertex_jit(path, i) == vertex_1 &&\n"             \
"            path_get_vertex_jit(path, i + 1) == vertex_2)\n"           \
"        {\n"                                                          \
"            return true;\n"                                            \
"        }\n"                                                          \
"    }\n"                                                              \
"    return false;\n"                                                   \
"}\n"

#define FIRST_MULTIPLE_PATHS_DEFN                                      \
"void first_multiple_paths_f(MultiplePaths *z, MultiplePaths *x,\n"     \
"    bool *_y)\n"                                                       \
"{\n"                                                                  \
"    (void) _y;\n"                                                      \
"    *z = *x;\n"                                                        \
"}\n"

#define SECOND_MULTIPLE_PATHS_DEFN                                     \
"void second_multiple_paths_f(MultiplePaths *z, bool *_x,\n"           \
"    MultiplePaths *y)\n"                                               \
"{\n"                                                                  \
"    (void) _x;\n"                                                      \
"    *z = *y;\n"                                                        \
"}\n"

#define COMBINE_MULTIPLE_PATHS_DEFN                                    \
"void combine_multiple_paths_f(MultiplePaths *z,\n"                    \
"    const MultiplePaths *x, const MultiplePaths *y)\n"                 \
"{\n"                                                                  \
"    MultiplePaths x_copy = *x;\n"                                      \
"    MultiplePaths y_copy = *y;\n"                                      \
"    size_t path_count = x_copy.path_count + y_copy.path_count;\n"      \
"    if (!multiple_paths_prepare_jit(z, path_count)) return;\n"         \
"    for (size_t j = 0; j < x_copy.path_count; j++)\n"                  \
"    {\n"                                                              \
"        multiple_paths_append_unchecked_jit(z,\n"                     \
"            multiple_paths_nth_const_jit(&x_copy, j));\n"             \
"    }\n"                                                              \
"    for (size_t j = 0; j < y_copy.path_count; j++)\n"                  \
"    {\n"                                                              \
"        multiple_paths_append_unchecked_jit(z,\n"                     \
"            multiple_paths_nth_const_jit(&y_copy, j));\n"             \
"    }\n"                                                              \
"}\n"

#define EXTEND_MULTIPLE_PATHS_DEFN                                     \
"void extend_multiple_paths_f(MultiplePaths *z,\n"                     \
"    const MultiplePaths *x, GrB_Index _row, GrB_Index col,\n"          \
"    const void *_y)\n"                                                 \
"{\n"                                                                  \
"    MultiplePaths src = *x;\n"                                         \
"    (void) _row;\n"                                                    \
"    (void) _y;\n"                                                      \
"    if (!multiple_paths_prepare_jit(z, src.path_count)) return;\n"     \
"    for (size_t i = 0; i < src.path_count; i++)\n"                     \
"    {\n"                                                              \
"        Path path = *multiple_paths_nth_const_jit(&src, i);\n"        \
"        if (!all_paths_can_extend_jit(&path)) continue;\n"            \
"        path_extend_jit(&path, (Vertex) col);\n"                      \
"        if (!path_is_empty_jit(&path))\n"                              \
"        {\n"                                                          \
"            multiple_paths_append_unchecked_jit(z, &path);\n"         \
"        }\n"                                                          \
"    }\n"                                                              \
"}\n"

#define EXTEND_MULTIPLE_SIMPLE_DEFN                                    \
"void extend_multiple_simple_f(MultiplePaths *z,\n"                    \
"    const MultiplePaths *x, GrB_Index _row, GrB_Index col,\n"          \
"    const void *_y)\n"                                                 \
"{\n"                                                                  \
"    MultiplePaths src = *x;\n"                                         \
"    (void) _row;\n"                                                    \
"    (void) _y;\n"                                                      \
"    if (!multiple_paths_prepare_jit(z, src.path_count)) return;\n"     \
"    for (size_t i = 0; i < src.path_count; i++)\n"                     \
"    {\n"                                                              \
"        Path path = *multiple_paths_nth_const_jit(&src, i);\n"        \
"        if (path_extending_will_add_repeated_vertex_jit(\n"            \
"            &path, (Vertex) col))\n"                                   \
"        {\n"                                                          \
"            continue;\n"                                               \
"        }\n"                                                          \
"        path_extend_jit(&path, (Vertex) col);\n"                      \
"        if (!path_is_empty_jit(&path))\n"                              \
"        {\n"                                                          \
"            multiple_paths_append_unchecked_jit(z, &path);\n"         \
"        }\n"                                                          \
"    }\n"                                                              \
"}\n"

#define EXTEND_MULTIPLE_TRAILS_DEFN                                    \
"void extend_multiple_trails_f(MultiplePaths *z,\n"                    \
"    const MultiplePaths *x, GrB_Index _row, GrB_Index col,\n"          \
"    const void *y)\n"                                                  \
"{\n"                                                                  \
"    MultiplePaths src = *x;\n"                                         \
"    (void) _row;\n"                                                    \
"    (void) y;\n"                                                       \
"    if (!multiple_paths_prepare_jit(z, src.path_count)) return;\n"     \
"    for (size_t i = 0; i < src.path_count; i++)\n"                     \
"    {\n"                                                              \
"        Path path = *multiple_paths_nth_const_jit(&src, i);\n"        \
"        if (path_extending_will_add_repeated_edge_jit(&path,\n"        \
"            (Vertex) col))\n"                                          \
"        {\n"                                                          \
"            continue;\n"                                               \
"        }\n"                                                          \
"        path_extend_jit(&path, (Vertex) col);\n"                      \
"        if (!path_is_empty_jit(&path))\n"                              \
"        {\n"                                                          \
"            multiple_paths_append_unchecked_jit(z, &path);\n"         \
"        }\n"                                                          \
"    }\n"                                                              \
"}\n"
//

// Helpers for Path with inline vertices + heap/temp overflow vertices.
static size_t path_extra_len (const Path *path)
{
    if (path->vertex_count > QUICK_PATH_LENGTH)
    {
        return path->vertex_count - QUICK_PATH_LENGTH ;
    }

    return 0 ;
}

static Vertex path_get_vertex (const Path *path, size_t i)
{
    assert (path != NULL) ;
    assert (i < path->vertex_count) ;

    if (i < QUICK_PATH_LENGTH)
    {
        return path->vertices[i] ;
    }

    assert (path->extra != NULL) ;
    return path->extra->vertices[i - QUICK_PATH_LENGTH] ;
}

static bool path_is_empty (const Path *path)
{
    return path->vertex_count == 0 ;
}

static Vertex path_start_vertex (const Path *path)
{
    return path_get_vertex (path, 0) ;
}

static Vertex path_last_vertex (const Path *path)
{
    return path_get_vertex (path, path->vertex_count - 1) ;
}

static bool path_is_closed_cycle (const Path *path)
{
    return path->vertex_count > 1 &&
        path_start_vertex (path) == path_last_vertex (path) ;
}

static bool all_paths_can_extend (const Path *path)
{
    return RPQ_MAX_PATH_LENGTH == 0 ||
        path->vertex_count < RPQ_MAX_PATH_LENGTH ;
}

static PathExtra *path_extra_temp_alloc_copy_plus_one
(
    const Path *src,
    Vertex vertex
)
{
    size_t old_extra_len ;
    size_t new_extra_len ;
    size_t nbytes ;
    PathExtra *extra ;

    old_extra_len = path_extra_len (src) ;
    new_extra_len = old_extra_len + 1 ;
    nbytes = sizeof (PathExtra) + new_extra_len * sizeof (Vertex) ;

    extra = (PathExtra *) temp_calloc_bytes (nbytes) ;
    if (extra == NULL)
    {
        return NULL ;
    }

    extra->len = new_extra_len ;

    if (old_extra_len > 0)
    {
        assert (src->extra != NULL) ;
        memcpy (extra->vertices, src->extra->vertices,
            old_extra_len * sizeof (Vertex)) ;
    }

    extra->vertices[old_extra_len] = vertex ;
    return extra ;
}

static void path_extend (Path *path, Vertex vertex)
{
    PathExtra *new_extra ;

    if (path->vertex_count == 0)
    {
        return ;
    }

    if (path->vertex_count < QUICK_PATH_LENGTH)
    {
        path->vertices[path->vertex_count] = vertex ;
        path->vertex_count++ ;
        return ;
    }

    new_extra = path_extra_temp_alloc_copy_plus_one (path, vertex) ;
    if (new_extra == NULL)
    {
        path->vertex_count = 0 ;
        path->extra = NULL ;
        return ;
    }

    path->extra = new_extra ;
    path->vertex_count++ ;
}

static int path_clone_heap (Path *dst, const Path *src, char *msg)
{
    size_t extra_len ;
    size_t nbytes ;
    int info ;

    memset (dst, 0, sizeof (*dst)) ;
    dst->vertex_count = src->vertex_count ;
    memcpy (dst->vertices, src->vertices, sizeof (dst->vertices)) ;

    extra_len = path_extra_len (src) ;
    if (extra_len == 0)
    {
        dst->extra = NULL ;
        return GrB_SUCCESS ;
    }

    nbytes = sizeof (PathExtra) + extra_len * sizeof (Vertex) ;
    info = LAGraph_Malloc ((void **) &dst->extra, nbytes, sizeof (char), msg) ;
    if (info != GrB_SUCCESS)
    {
        return info ;
    }

    dst->extra->len = extra_len ;
    memcpy (dst->extra->vertices, src->extra->vertices,
        extra_len * sizeof (Vertex)) ;
    return GrB_SUCCESS ;
}

static void path_destroy_heap (Path *path)
{
    if (path == NULL)
    {
        return ;
    }

    LAGraph_Free ((void **) &path->extra, NULL) ;
    path->vertex_count = 0 ;
}

static void free_all (Path **paths, size_t path_count)
{
    if (paths == NULL || *paths == NULL)
    {
        return ;
    }

    for (size_t i = 0 ; i < path_count ; i++)
    {
        path_destroy_heap (&((*paths)[i])) ;
    }

    LAGraph_Free ((void **) paths, NULL) ;
}
//

void Path_print (const Path *x)
{
    if (x->vertex_count == 0)
    {
        printf ("empty path \n") ;
        return ;
    }

    for (size_t i = 0 ; i < x->vertex_count ; i++)
    {
        // Increase the vertex by 1 since usually user expects the same
        // numbering as in the input determined by MTX file in which the
        // entries are enumerated starting from 1.
        printf ("(%llu)", (unsigned long long) (path_get_vertex (x, i) + 1)) ;

        if (i != x->vertex_count - 1)
        {
            printf ("-") ;
        }
    }

    printf ("\n") ;
}

// Block with helper function for handling cases
// when we have more paths than QUICK_PATH_COUNT
static size_t multiple_paths_capacity (const MultiplePaths *x)
{
    size_t cap = QUICK_PATH_COUNT ;

    if (x->extra != NULL)
    {
        cap += x->extra->count ;
    }

    return cap ;
}

static const Path *multiple_paths_nth_const
(
    const MultiplePaths *x,
    size_t j
)
{
    assert (j < x->path_count) ;

    if (j < QUICK_PATH_COUNT)
    {
        return &x->paths[j] ;
    }

    assert (x->extra != NULL) ;
    return &x->extra->paths[j - QUICK_PATH_COUNT] ;
}

static Path *multiple_paths_nth_mut
(
    MultiplePaths *x,
    size_t j
)
{
    assert (j < multiple_paths_capacity (x)) ;

    if (j < QUICK_PATH_COUNT)
    {
        return &x->paths[j] ;
    }

    assert (x->extra != NULL) ;
    return &x->extra->paths[j - QUICK_PATH_COUNT] ;
}

static MultiplePathsExtra *multiple_paths_extra_temp_alloc (size_t extra_count)
{
    size_t nbytes ;
    MultiplePathsExtra *extra ;

    if (extra_count == 0)
    {
        return NULL ;
    }

    nbytes = sizeof (MultiplePathsExtra) + extra_count * sizeof (Path) ;
    extra = (MultiplePathsExtra *) temp_calloc_bytes (nbytes) ;
    if (extra == NULL)
    {
        return NULL ;
    }

    extra->count = extra_count ;
    return extra ;
}

static bool multiple_paths_prepare (MultiplePaths *x, size_t capacity)
{
    memset (x, 0, sizeof (*x)) ;

    // Preserve the per-point path cap. If one frontier cell
    // would contain too many paths, GraphBLAS UDF cannot return an error
    // directly, so we mark a global flag and check it after GrB_mxm/GrB_apply.
    if (capacity > PATHS_PER_POINT_LIMIT)
    {
        path_limit_mark_exceeded () ;
        return false ;
    }

    if (capacity > QUICK_PATH_COUNT)
    {
        x->extra = multiple_paths_extra_temp_alloc (capacity - QUICK_PATH_COUNT) ;
        if (x->extra == NULL)
        {
            x->path_count = 0 ;
            return false ;
        }
    }

    return true ;
}

static bool multiple_paths_append (MultiplePaths *x, const Path *path)
{
    // handle per-point cap
    if (x->path_count >= PATHS_PER_POINT_LIMIT)
    {
        path_limit_mark_exceeded () ;
        return false ;
    }

    if (x->path_count >= multiple_paths_capacity (x))
    {
        return false ;
    }

    *multiple_paths_nth_mut (x, x->path_count) = *path ;
    x->path_count++ ;
    return true ;
}

static void multiple_paths_append_unchecked (MultiplePaths *x, const Path *path)
{
    *multiple_paths_nth_mut (x, x->path_count) = *path ;
    x->path_count++ ;
}
//

// All functions below reworked.
// Due to graphblas api, we must handle z param like it's empty.
// It should just store result (z = f(x)).
// So we can't use it for any checks, or use its fields for something.
void first_multiple_paths_f(MultiplePaths *z, MultiplePaths *x, bool *_y)
{
    (void) _y ;
    *z = *x ;
}

void second_multiple_paths_f(MultiplePaths *z, bool *_x, MultiplePaths *y)
{
    (void) _x ;
    *z = *y ;
}

void combine_multiple_paths_f(MultiplePaths *z, const MultiplePaths *x, const MultiplePaths *y)
{
    MultiplePaths x_copy = *x ;
    MultiplePaths y_copy = *y ;
    size_t path_count = x_copy.path_count + y_copy.path_count ;

    if (!multiple_paths_prepare (z, path_count))
    {
        return ;
    }

    for (size_t j = 0 ; j < x_copy.path_count ; j++)
    {
        multiple_paths_append_unchecked (z, multiple_paths_nth_const (&x_copy, j)) ;
    }

    for (size_t j = 0 ; j < y_copy.path_count ; j++)
    {
        multiple_paths_append_unchecked (z, multiple_paths_nth_const (&y_copy, j)) ;
    }
}

//
// ALL PATHS.
//

// NB: Using this semantic without a length limit makes the code behave like a
// procedure for searching all paths satisfying the constraints.
// It means it may not finish if there is loops.

void extend_multiple_paths_f(MultiplePaths *z, const MultiplePaths *x, GrB_Index _row, GrB_Index col, const void *_y)
{
    MultiplePaths src = *x ;

    (void) _row ;
    (void) _y ;

    if (!multiple_paths_prepare (z, src.path_count))
    {
        return ;
    }

    for (size_t i = 0 ; i < src.path_count ; i++)
    {
        Path path = *multiple_paths_nth_const (&src, i) ;
        if (!all_paths_can_extend (&path))
        {
            continue ;
        }

        path_extend (&path, (Vertex) col) ;

        if (!path_is_empty (&path))
        {
            multiple_paths_append_unchecked (z, &path) ;
        }
    }
}

//
// ALL SIMPLE
//


// remove cycle support from all simple
static inline bool path_extending_will_add_repeated_vertex(const Path *path, Vertex vertex)
{
    if (path->vertex_count == 0)
    {
        return false ;
    }

    if (path_is_closed_cycle (path))
    {
        return true ;
    }

    for (size_t i = 0 ; i < path->vertex_count ; i++)
    {
        if (path_get_vertex (path, i) == vertex)
        {
            return true ;
        }
    }

    return false ;
}

void extend_multiple_simple_f(MultiplePaths *z, const MultiplePaths *x, GrB_Index _row, GrB_Index col, const void *_y)
{
    MultiplePaths src = *x ;

    (void) _row ;
    (void) _y ;

    if (!multiple_paths_prepare (z, src.path_count))
    {
        return ;
    }

    for (size_t i = 0 ; i < src.path_count ; i++)
    {
        Path path = *multiple_paths_nth_const (&src, i) ;

        if (path_extending_will_add_repeated_vertex (&path,
            (Vertex) col))
        {
            continue ;
        }

        path_extend (&path, (Vertex) col) ;
        if (!path_is_empty (&path))
        {
            multiple_paths_append_unchecked (z, &path) ;
        }
    }
}

//
// ALL TRAILS
//

static inline bool path_extending_will_add_repeated_edge(const Path *path, Vertex vertex_2)
{
    if (path->vertex_count == 0)
    {
        return false ;
    }

    // We identify edges as pairs of vertices.
    Vertex vertex_1 = path_last_vertex (path) ;

    for (size_t i = 0 ; i + 1 < path->vertex_count ; i++)
    {
        if (path_get_vertex (path, i) == vertex_1 &&
            path_get_vertex (path, i + 1) == vertex_2)
        {
            return true ;
        }
    }

    return false ;
}

void extend_multiple_trails_f(MultiplePaths *z, const MultiplePaths *x, GrB_Index _row, GrB_Index col, const void *y)
{
    MultiplePaths src = *x ;

    (void) _row ;
    (void) y ;

    if (!multiple_paths_prepare (z, src.path_count))
    {
        return ;
    }

    for (size_t i = 0 ; i < src.path_count ; i++)
    {
        Path path = *multiple_paths_nth_const (&src, i) ;

        if (path_extending_will_add_repeated_edge (&path, (Vertex) col))
        {
            continue ;
        }

        path_extend (&path, (Vertex) col) ;
        if (!path_is_empty (&path))
        {
            multiple_paths_append_unchecked (z, &path) ;
        }
    }
}

// Result array can grow past PATH_LIMIT if the caller passes a larger limit.
static int ensure_result_capacity
(
    Path **paths,
    size_t *capacity,
    size_t need,
    char *msg
)
{
    Path *new_paths ;
    size_t new_capacity ;
    int info ;

    if (need <= *capacity)
    {
        return GrB_SUCCESS ;
    }

    new_capacity = *capacity ;
    if (new_capacity == 0)
    {
        new_capacity = PATH_LIMIT ;
        if (new_capacity == 0)
        {
            new_capacity = 1 ;
        }
    }

    while (new_capacity < need)
    {
        if (new_capacity > SIZE_MAX / 2)
        {
            new_capacity = need ;
            break ;
        }

        new_capacity = new_capacity * 2 ;
    }

    info = LAGraph_Malloc ((void **) &new_paths, new_capacity,
        sizeof (Path), msg) ;
    if (info != GrB_SUCCESS)
    {
        return info ;
    }

    memset (new_paths, 0, new_capacity * sizeof (Path)) ;

    if (*paths != NULL)
    {
        memcpy (new_paths, *paths, (*capacity) * sizeof (Path)) ;
        LAGraph_Free ((void **) paths, NULL) ;
    }

    *paths = new_paths ;
    *capacity = new_capacity ;
    return GrB_SUCCESS ;
}

static int final_state_was_visited
(
    GrB_Matrix visited,
    const GrB_Index *QF,
    size_t nqf,
    GrB_Index vertex,
    bool *seen
)
{
    GrB_Info info ;
    bool value ;

    *seen = false ;

    for (size_t i = 0 ; i < nqf ; i++)
    {
        value = false ;
        info = GrB_Matrix_extractElement_BOOL (&value, visited, QF [i],
            vertex) ;
        if (info == GrB_NO_VALUE)
        {
            continue ;
        }
        if (info != GrB_SUCCESS)
        {
            return info ;
        }
        if (value)
        {
            *seen = true ;
            return GrB_SUCCESS ;
        }
    }

    return GrB_SUCCESS ;
}


#undef LG_FREE_WORK
#undef LG_FREE_ALL
#define LG_FREE_WORK                            \
{                                               \
    GrB_free (&frontier) ;                      \
    GrB_free (&next_frontier) ;                 \
    GrB_free (&symbol_frontier) ;               \
    GrB_free (&visited) ;                       \
    LAGraph_Free ((void **) &A, NULL) ;         \
    LAGraph_Free ((void **) &AT, NULL) ;        \
    LAGraph_Free ((void **) &B, NULL) ;         \
    LAGraph_Free ((void **) &BT, NULL) ;        \
    LAGraph_Free ((void **) &X, NULL) ;         \
    LAGraph_Free ((void **) &I, NULL) ;         \
    LAGraph_Free ((void **) &J, NULL) ;         \
    temp_arena_destroy () ;                     \
}

#define LG_FREE_ALL                             \
{                                               \
    LG_FREE_WORK ;                              \
    free_all (paths, *path_count) ;             \
    if (paths != NULL) *paths = NULL ;          \
    if (path_count != NULL) *path_count = 0 ;   \
}

static int LAGraph_2Rpq
(
    // output:
    Path **paths,               // simple paths from one of the starting
                                // nodes satisfying regular constraints
    size_t *path_count,         // resulting path count
    // input:
    LAGraph_Graph *R,           // input non-deterministic finite automaton
                                // adjacency matrix decomposition
    bool *inverse_labels,       // inversed labels
    size_t nl,                  // total label count, # of matrices graph and
                                // NFA adjacency matrix decomposition
    const GrB_Index *QS,        // starting states in NFA
    size_t nqs,                 // number of starting states in NFA
    const GrB_Index *QF,        // final states in NFA
    size_t nqf,                 // number of final states in NFA
    LAGraph_Graph *G,           // input graph adjacency matrix decomposition
    const GrB_Index *S,         // source vertices to start searching paths
    size_t ns,                  // number of source vertices
    bool inverse,               // inverse the whole query
    bool ignore_visited,        // use mask to avoid processing the same (q, v)
    uint64_t limit,             // maximum path count
    char *msg,                  // LAGraph output message
    GrB_IndexUnaryOp op         // index unary op for a specific semantic
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;

    GrB_Matrix frontier = NULL ;         // traversal frontier representing
                                         // correspondence between NFA states
                                         // and graph vertices
    GrB_Matrix symbol_frontier = NULL ;  // part of the new frontier for the
                                         // specific label
    GrB_Matrix next_frontier = NULL ;    // frontier value on the next
                                         // traversal step
    GrB_Matrix visited = NULL ;          // visited pairs (state, vertex)

    GrB_Index ng = 0 ;                   // # nodes in the graph
    GrB_Index nr = 0 ;                   // # states in the NFA

    GrB_Index rows = 0 ;                 // utility matrix row count
    GrB_Index cols = 0 ;                 // utility matrix column count

    // TODO: This names might be too short.
    GrB_Semiring sr1 = first_combine_multiple_paths ;
    GrB_Semiring sr2 = second_combine_multiple_paths ;
    GrB_BinaryOp acc = combine_multiple_paths_op ;

    GrB_Matrix *A = NULL ;
    GrB_Matrix *AT = NULL ;
    GrB_Matrix *B = NULL ;
    GrB_Matrix *BT = NULL ;

    MultiplePaths *X = NULL ;
    GrB_Index *I = NULL ;
    GrB_Index *J = NULL ;
    size_t result_capacity = 0 ;

    if (paths == NULL || path_count == NULL || G == NULL || R == NULL ||
        S == NULL || op == NULL)
    {
        return GrB_NULL_POINTER ;
    }

    (*paths) = NULL ;
    (*path_count) = 0 ;

    path_limit_reset () ;

    LG_ASSERT_MSG (!ignore_visited || ns == 1, GrB_INVALID_VALUE,
        "AllShortestPaths requires exactly one source vertex") ;

    // init arenas for pathExtra
    LG_TRY (temp_arena_init (msg)) ;

    for (size_t i = 0 ; i < nl ; i++)
    {
        if (G[i] == NULL) continue ;
        LG_TRY (LAGraph_CheckGraph (G[i], msg)) ;
    }

    for (size_t i = 0 ; i < nl ; i++)
    {
        if (R[i] == NULL) continue ;
        LG_TRY (LAGraph_CheckGraph (R[i], msg)) ;
    }

    LG_TRY (LAGraph_Malloc ((void **) &A, nl, sizeof (GrB_Matrix), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &AT, nl, sizeof (GrB_Matrix), msg)) ;

    for (size_t i = 0 ; i < nl ; i++)
    {
        if (G[i] == NULL)
        {
            A[i] = NULL ;
            AT[i] = NULL ;
            continue ;
        }

        A[i] = G[i]->A ;
        if (G[i]->kind == LAGraph_ADJACENCY_UNDIRECTED ||
            G[i]->is_symmetric_structure == LAGraph_TRUE)
        {
            AT[i] = A[i] ;
        }
        else
        {
            // AT[i] could be NULL and the matrix will be transposed by a
            // descriptor
            AT[i] = G[i]->AT ;
        }
    }

    LG_TRY (LAGraph_Malloc ((void **) &B, nl, sizeof (GrB_Matrix), msg)) ;
    LG_TRY (LAGraph_Malloc ((void **) &BT, nl, sizeof (GrB_Matrix), msg)) ;

    for (size_t i = 0 ; i < nl ; i++)
    {
        BT[i] = NULL ;

        if (R[i] == NULL)
        {
            B[i] = NULL ;
            BT[i] = NULL ;
            continue ;
        }

        B[i] = R[i]->A ;
        if (R[i]->is_symmetric_structure == LAGraph_TRUE)
        {
            BT[i] = B[i] ;
        }
        else
        {
            // BT[i] could be NULL and the matrix will be transposed by a
            // descriptor
            BT[i] = R[i]->AT ;
        }
    }

    for (size_t i = 0 ; i < nl ; i++)
    {
        if (A[i] == NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&ng, A[i])) ;
        break ;
    }

    for (size_t i = 0 ; i < nl ; i++)
    {
        if (B[i] == NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&nr, B[i])) ;
        break ;
    }

    // Check all the matrices in graph adjacency matrix decomposition are
    // square and of the same dimensions
    for (size_t i = 0 ; i < nl ; i++)
    {
        if (A[i] == NULL) continue ;

        GRB_TRY (GrB_Matrix_nrows (&rows, A[i])) ;
        GRB_TRY (GrB_Matrix_ncols (&cols, A[i])) ;

        LG_ASSERT_MSG (rows == ng && cols == ng, LAGRAPH_NOT_CACHED,
            "all the matrices in the graph adjacency matrix decomposition "
            "should have the same dimensions and be square") ;
    }

    // Check all the matrices in NFA adjacency matrix decomposition are
    // square and of the same dimensions
    for (size_t i = 0 ; i < nl ; i++)
    {
        if (B[i] == NULL) continue ;

        GrB_Index rrows = 0 ;
        GrB_Index rcols = 0 ;

        GRB_TRY (GrB_Matrix_nrows (&rrows, B[i])) ;
        GRB_TRY (GrB_Matrix_ncols (&rcols, B[i])) ;

        LG_ASSERT_MSG (rrows == nr && rcols == nr, LAGRAPH_NOT_CACHED,
            "all the matrices in the NFA adjacency matrix decomposition "
            "should have the same dimensions and be square") ;
    }

    // Check source nodes in the graph
    for (size_t i = 0 ; i < ns ; i++)
    {
        GrB_Index s = S [i] ;
        LG_ASSERT_MSG (s < ng, GrB_INVALID_INDEX, "invalid graph source node") ;
    }

    // Check starting states of the NFA
    for (size_t i = 0 ; i < nqs ; i++)
    {
        GrB_Index qs = QS [i] ;
        LG_ASSERT_MSG (qs < nr, GrB_INVALID_INDEX,
            "invalid NFA starting state") ;
    }

    // Check final states of the NFA
    for (size_t i = 0 ; i < nqf ; i++)
    {
        GrB_Index qf = QF [i] ;
        LG_ASSERT_MSG (qf < nr, GrB_INVALID_INDEX, "invalid NFA final state") ;
    }

    // -------------------------------------------------------------------------
    // initialization
    // -------------------------------------------------------------------------

    GRB_TRY (LAGraph_Calloc ((void **) paths, PATH_LIMIT,
        sizeof (Path), msg)) ;
    result_capacity = PATH_LIMIT ;

    GRB_TRY (GrB_Matrix_new (&next_frontier, multiple_paths, nr, ng)) ;

    if (ignore_visited)
    {
        GRB_TRY (GrB_Matrix_new (&visited, GrB_BOOL, nr, ng)) ;
    }

    // Initialize frontier with the source nodes

    for (size_t i = 0 ; i < ns ; i++)
    {
        GrB_Index s = S[i] ;
        MultiplePaths value = {
            .paths = {
                {
                    .vertices = { s },
                    .vertex_count = 1,
                    .extra = NULL
                }
            },
            .path_count = 1,
            .extra = NULL
        };

        for (size_t j = 0 ; j < nqs ; j++)
        {
            GrB_Index qs = QS[j] ;

            GRB_TRY (GrB_Matrix_setElement_UDT (next_frontier, &value, qs, s)) ;
        }
    }

    // Initialize a few utility matrices
    GRB_TRY (GrB_Matrix_new (&frontier, multiple_paths, nr, ng)) ;
    GRB_TRY (GrB_Matrix_new (&symbol_frontier, multiple_paths, nr, ng)) ;

    // Main loop
    while (true)
    {
        GrB_Index nvals = 0 ;
        bool had_non_empty_path = false ;

        GRB_TRY (GrB_Matrix_nvals (&nvals, next_frontier)) ;

        LG_TRY (LAGraph_Calloc ((void **) &X, nvals, sizeof (MultiplePaths), msg)) ;
        LG_TRY (LAGraph_Calloc ((void **) &I, nvals, sizeof (GrB_Index), msg)) ;
        if (ignore_visited)
        {
            LG_TRY (LAGraph_Calloc ((void **) &J, nvals,
                sizeof (GrB_Index), msg)) ;
        }

        // TODO: Change to a generic call.
        GRB_TRY (GrB_Matrix_extractTuples_UDT (I,
            ignore_visited ? J : GrB_NULL, (void**) X, &nvals,
            next_frontier)) ;

        for (size_t i = 0 ; i < nvals ; i++)
        {
            for (size_t j = 0 ; j < X[i].path_count ; j++)
            {
                const Path *path = multiple_paths_nth_const (&X[i], j) ;

                // Required beacause we need to handle not only quick paths
                if (!path_is_empty (path))
                {
                    had_non_empty_path = true;
                    break;
                }
            }

            bool final = false ;
            for (size_t j = 0 ; j < nqf ; j++)
            {
                if (I[i] == QF[j])
                {
                    final = true ;
                    break ;
                }
            }

            if (!final)
            {
                continue ;
            }

            if (ignore_visited)
            {
                bool seen_final = false ;

                GRB_TRY (final_state_was_visited (visited, QF, nqf, J [i],
                    &seen_final)) ;
                if (seen_final)
                {
                    continue ;
                }
            }

            if ((*path_count) >= limit)
            {
                continue ;
            }

            uint64_t remaining_limit = limit - (uint64_t) (*path_count) ;
            size_t paths_to_reserve = X[i].path_count ;
            if (remaining_limit < (uint64_t) paths_to_reserve)
            {
                paths_to_reserve = (size_t) remaining_limit ;
            }
            size_t result_need = (*path_count) + paths_to_reserve ;
            if (result_need > result_capacity)
            {
                LG_TRY (ensure_result_capacity (paths, &result_capacity,
                    result_need, msg)) ;
            }

            for (size_t j = 0 ; j < X[i].path_count && (*path_count) < limit ; j++)
            {
                // Deep copy of result
                const Path *path = multiple_paths_nth_const (&X[i], j) ;
                if (!path_is_empty (path))
                {
                    LG_TRY (path_clone_heap (&((*paths)[(*path_count)]), path,
                        msg)) ;
                    (*path_count)++ ;
                }
            }
        }

        if (ignore_visited)
        {
            GRB_TRY (GrB_assign (visited, next_frontier, GrB_NULL, true,
                GrB_ALL, nr, GrB_ALL, ng, GrB_DESC_S)) ;
        }

        LAGraph_Free ((void **) &X, NULL) ;
        LAGraph_Free ((void **) &I, NULL) ;
        LAGraph_Free ((void **) &J, NULL) ;

        if (!had_non_empty_path || (*path_count) == limit)
        {
            break;
        }

        GrB_Matrix old_frontier = frontier ;
        frontier = next_frontier ;
        next_frontier = old_frontier ;

        GRB_TRY (GrB_Matrix_clear(next_frontier)) ;

        // Obtain a new relation between the NFA states and the graph nodes
        for (size_t i = 0 ; i < nl ; i++)
        {
            if (A[i] == NULL || B[i] == NULL) continue ;

            // Traverse the NFA
            // Try to use a provided transposed matrix or use the descriptor
            if (!inverse) {
                if (BT[i] != NULL)
                {
                    GRB_TRY (GrB_mxm (symbol_frontier, GrB_NULL, GrB_NULL,
                        sr2, BT[i], frontier, GrB_DESC_R)) ;
                }
                else
                {
                    GRB_TRY (GrB_mxm (symbol_frontier, GrB_NULL, GrB_NULL,
                        sr2, B[i], frontier, GrB_DESC_RT0)) ;
                }
            } else {
                GRB_TRY (GrB_mxm (symbol_frontier, GrB_NULL, GrB_NULL, sr2, B[i], frontier, GrB_DESC_R )) ;
            }

            // Skip the iteration if symbol_frontier is already empty.
            GrB_Index symbol_nvals = 0 ;
            GRB_TRY (GrB_Matrix_nvals (&symbol_nvals, symbol_frontier)) ;
            if (symbol_nvals == 0) continue ;

            GrB_Matrix mask = ignore_visited ? visited : GrB_NULL ;
            GrB_Descriptor desc_forward = ignore_visited ? GrB_DESC_SC :
                GrB_NULL ;
            GrB_Descriptor desc_backward = ignore_visited ? GrB_DESC_SCT1 :
                GrB_DESC_T1 ;

            // Traverse the graph
            if (!inverse_labels[i]) {
                if (!inverse) {
                    GRB_TRY (GrB_mxm (next_frontier, mask, acc, sr1, symbol_frontier, A[i], desc_forward)) ;
                } else if (AT[i]) {
                    GRB_TRY (GrB_mxm (next_frontier, mask, acc, sr1, symbol_frontier, AT[i], desc_forward)) ;
                } else {
                    GRB_TRY (GrB_mxm (next_frontier, mask, acc, sr1, symbol_frontier, A[i], desc_backward)) ;
                }
            } else {
                if (!inverse && AT[i]) {
                    GRB_TRY (GrB_mxm (next_frontier, mask, acc, sr1, symbol_frontier, AT[i], desc_forward)) ;
                } else if (!inverse) {
                    GRB_TRY (GrB_mxm (next_frontier, mask, acc, sr1, symbol_frontier, A[i], desc_backward)) ;
                } else {
                    GRB_TRY (GrB_mxm (next_frontier, mask, acc, sr1, symbol_frontier, A[i], desc_forward)) ;
                }
            }

            // Handle oom
            LG_ASSERT_MSGF (!temp_alloc_failed (), GrB_OUT_OF_MEMORY,
                "out of memory in temporary RPQ allocator: used=%zu capacity=%zu",
                temp_arena_used (), temp_arena_capacity) ;
            if (path_limit_failed ())
            {
                LG_ERROR_MSG ("LAGraph failure (file %s, line %d): %s",
                    __FILE__, __LINE__, "path limit per point exceeded") ;
                LG_FREE_WORK ;
                return GrB_OUT_OF_MEMORY ;
            }
        }

        GRB_TRY (GrB_apply (next_frontier, GrB_NULL, GrB_NULL, op, next_frontier, false, GrB_NULL)) ;

        // Handle oom
        LG_ASSERT_MSGF (!temp_alloc_failed (), GrB_OUT_OF_MEMORY,
            "out of memory in temporary RPQ allocator: used=%zu capacity=%zu",
            temp_arena_used (), temp_arena_capacity) ;
        if (path_limit_failed ())
        {
            LG_ERROR_MSG ("LAGraph failure (file %s, line %d): %s",
                __FILE__, __LINE__, "path limit per point exceeded") ;
            LG_FREE_WORK ;
            return GrB_OUT_OF_MEMORY ;
        }
    }

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}


int LAGraph_2Rpq_AllSimple      // All simple paths satisfying regular
                                // expression. Simple paths are paths without
                                // repeated vertices.
(
    // output:
    Path **paths,               // simple paths from one of the starting
                                // nodes satisfying regular constraints
    size_t *path_count,         // resulting path count
    // input:
    LAGraph_Graph *R,           // input non-deterministic finite automaton
                                // adjacency matrix decomposition
    bool *inverse_labels,       // inversed labels
    size_t nl,                  // total label count, # of matrices graph and
                                // NFA adjacency matrix decomposition
    const GrB_Index *QS,        // starting states in NFA
    size_t nqs,                 // number of starting states in NFA
    const GrB_Index *QF,        // final states in NFA
    size_t nqf,                 // number of final states in NFA
    LAGraph_Graph *G,           // input graph adjacency matrix decomposition
    const GrB_Index *S,         // source vertices to start searching paths
    size_t ns,                  // number of source vertices
    bool inverse,               // inverse the whole query
    char *msg                   // LAGraph output message
)
{
    return LAGraph_2Rpq(paths, path_count, R, inverse_labels, nl, QS, nqs, QF, nqf, G, S, ns, inverse, false, ULLONG_MAX, msg, extend_multiple_simple) ;
}

LAGRAPHX_PUBLIC
int LAGraph_2Rpq_AllTrails      // All trails satisfying regular expression.
                                // Trails are paths without repeated edges.
(
    // output:
    Path **paths,                // trails from one of the starting nodes
                                // satisfying regular constraints
    size_t *path_count,         // resulting path count
    // input:
    LAGraph_Graph *R,           // input non-deterministic finite automaton
                                // adjacency matrix decomposition
    bool *inverse_labels,       // inversed labels
    size_t nl,                  // total label count, # of matrices graph and
                                // NFA adjacency matrix decomposition
    const GrB_Index *QS,        // starting states in NFA
    size_t nqs,                 // number of starting states in NFA
    const GrB_Index *QF,        // final states in NFA
    size_t nqf,                 // number of final states in NFA
    LAGraph_Graph *G,           // input graph adjacency matrix decomposition
    const GrB_Index *S,         // source vertices to start searching paths
    size_t ns,                  // number of source vertices
    bool inverse,               // inverse the whole query
    char *msg                   // LAGraph output message
)
{
    return LAGraph_2Rpq(paths, path_count, R, inverse_labels, nl, QS, nqs, QF, nqf, G, S, ns, inverse, false, ULLONG_MAX, msg, extend_multiple_trails) ;
}

int LAGraph_2Rpq_AllPaths       // All paths satisfying regular expression
(
    // output:
    Path **paths,               // paths from one of the starting nodes
                                // satisfying regular constraints
    size_t *path_count,         // resulting path count
                                // input:
    LAGraph_Graph *R,           // input non-deterministic finite automaton
                                // adjacency matrix decomposition
    bool *inverse_labels,       // inversed labels
    size_t nl,                  // total label count, # of matrices graph and
                                // NFA adjacency matrix decomposition
    const GrB_Index *QS,        // starting states in NFA
    size_t nqs,                 // number of starting states in NFA
    const GrB_Index *QF,        // final states in NFA
    size_t nqf,                 // number of final states in NFA
    LAGraph_Graph *G,           // input graph adjacency matrix decomposition
    const GrB_Index *S,         // source vertices to start searching paths
    size_t ns,                  // number of source vertices
    bool inverse,               // inverse the whole query
    uint64_t limit,             // maximum path count
    char *msg                   // LAGraph output message
    )
{
        return LAGraph_2Rpq(paths, path_count, R, inverse_labels, nl, QS, nqs, QF, nqf, G, S, ns, inverse, false, limit, msg, extend_multiple_paths) ;
}

int LAGraph_2Rpq_AllShortestPaths       // All shortest paths satisfying regular expression
(
    // output:
    Path **paths,               // paths from one of the starting nodes
                                // satisfying regular constraints
    size_t *path_count,         // resulting path count
                                // input:
    LAGraph_Graph *R,           // input non-deterministic finite automaton
                                // adjacency matrix decomposition
    bool *inverse_labels,       // inversed labels
    size_t nl,                  // total label count, # of matrices graph and
                                // NFA adjacency matrix decomposition
    const GrB_Index *QS,        // starting states in NFA
    size_t nqs,                 // number of starting states in NFA
    const GrB_Index *QF,        // final states in NFA
    size_t nqf,                 // number of final states in NFA
    LAGraph_Graph *G,           // input graph adjacency matrix decomposition
    const GrB_Index *S,         // source vertices to start searching paths
    size_t ns,                  // number of source vertices
    bool inverse,               // inverse the whole query
    uint64_t limit,             // maximum path count
    char *msg                   // LAGraph output message
    )
{
        return LAGraph_2Rpq(paths, path_count, R, inverse_labels, nl, QS, nqs, QF, nqf, G, S, ns, inverse, true, limit, msg, extend_multiple_paths) ;
}

// Required because returned Path objects may own heap-allocated PathExtra.
int LAGraph_2Rpq_FreePaths
(
    Path **paths,
    size_t path_count,
    char *msg
)
{
    (void) msg ;

    free_all (paths, path_count) ;
    return GrB_SUCCESS ;
}

#undef LG_FREE_WORK
#undef LG_FREE_ALL
#define LG_FREE_WORK                            \
{                                               \
}
#define LG_FREE_ALL                         \
{                                           \
        LG_FREE_WORK ;                          \
}

int LAGraph_Rpq_initialize(char *msg)
{
        (void) msg ;

        memset (&multiple_paths_identity, 0, sizeof (multiple_paths_identity)) ;

        GRB_TRY (GxB_Type_new (&multiple_paths, sizeof (MultiplePaths), "MultiplePaths", MULTIPLE_PATHS_TYPE_DEFN)) ;

        GRB_TRY (GxB_BinaryOp_new (&combine_multiple_paths_op, (GxB_binary_function) &combine_multiple_paths_f, multiple_paths, multiple_paths, multiple_paths, "combine_multiple_paths_f", COMBINE_MULTIPLE_PATHS_DEFN)) ;
        GRB_TRY (GxB_BinaryOp_new (&first_multiple_paths, (GxB_binary_function) &first_multiple_paths_f, multiple_paths, multiple_paths, GrB_BOOL, "first_multiple_paths_f", FIRST_MULTIPLE_PATHS_DEFN)) ;
        GRB_TRY (GxB_BinaryOp_new (&second_multiple_paths, (GxB_binary_function) &second_multiple_paths_f, multiple_paths, GrB_BOOL, multiple_paths, "second_multiple_paths_f", SECOND_MULTIPLE_PATHS_DEFN)) ;

        GRB_TRY (GrB_Monoid_new (&combine_multiple_paths, combine_multiple_paths_op, (void*) &multiple_paths_identity)) ;
        GRB_TRY (GrB_Semiring_new (&first_combine_multiple_paths, combine_multiple_paths, first_multiple_paths)) ;
        GRB_TRY (GrB_Semiring_new (&second_combine_multiple_paths, combine_multiple_paths, second_multiple_paths)) ;

        GRB_TRY (GxB_IndexUnaryOp_new (&extend_multiple_paths, (GxB_index_unary_function) &extend_multiple_paths_f, multiple_paths, multiple_paths, GrB_BOOL, "extend_multiple_paths_f", EXTEND_MULTIPLE_PATHS_DEFN)) ;
        GRB_TRY (GxB_IndexUnaryOp_new (&extend_multiple_simple, (GxB_index_unary_function) &extend_multiple_simple_f, multiple_paths, multiple_paths, GrB_BOOL, "extend_multiple_simple_f", EXTEND_MULTIPLE_SIMPLE_DEFN)) ;
        GRB_TRY (GxB_IndexUnaryOp_new (&extend_multiple_trails, (GxB_index_unary_function) &extend_multiple_trails_f, multiple_paths, multiple_paths, GrB_BOOL, "extend_multiple_trails_f", EXTEND_MULTIPLE_TRAILS_DEFN)) ;

        return GrB_SUCCESS ;
}
