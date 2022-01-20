//------------------------------------------------------------------------------
// LAGraph_Type: return the type or name of type of a matrix, vector, or scalar
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// On input, "char *name" is a pointer to a pre-allocated array of size at
// least LAGRAPH_MAX_NAME_LEN.  On output, the array is filled with a string
// corresponding to the type of a GrB_Matrix, GrB_Vector, or GrB_Scalar.
// For built-in types, the strings are defined as:

//      "bool"      GrB_BOOL
//      "int8_t"    GrB_INT8
//      "int16_t"   GrB_INT16
//      "int32_t"   GrB_INT32
//      "int64_t"   GrB_INT64
//      "uint8_t"   GrB_UINT8
//      "uint16_t"  GrB_UINT16
//      "uint32_t"  GrB_UINT32
//      "uint64_t"  GrB_UINT64
//      "float"     GrB_FP32
//      "double"    GrB_FP64

// SuiteSparse:GraphBLAS adds two extended types:
//      "float complex"     GxB_FC32
//      "double complex"    GxB_FC64

// For user-defined types, if SuiteSparse:GraphBLAS is used, then GrB_Type_new
// can capture the type name, if called as follows, where the 2nd parameter has
// the form "sizeof (T)" for some C typedef type T.
//
//      typedef ... myctype ;
//      GrB_Type MyType ;
//      GrB_Type_new (&MyType, sizeof (myctype)) ;
//
// In this case, LAGraph_*TypeName returns the string "myctype".

// Currently, these methods require SuiteSparse:GraphBLAS.  Other GraphBLAS
// libraries will result in a return value of GrB_NOT_IMPLEMENTED, and the name
// is returned as an empty string.  The type cannot be queried using the v2.0 C
// API.  This will be resolved in a future C API spec.

#include "LG_internal.h"

//------------------------------------------------------------------------------
// LAGraph_MatrixTypeName: return the name of the GrB_Type of a GrB_Matrix
//------------------------------------------------------------------------------

int LAGraph_MatrixTypeName // returns 0 if successful, < 0 if failure
(
    char *name,             // name of the type of the matrix A
    GrB_Matrix A,           // matrix to query
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (name != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // determine the name of the type of the GrB_Matrix A
    //--------------------------------------------------------------------------

    #if LG_SUITESPARSE
    return (GxB_Matrix_type_name (name, A)) ;
    #else
    name [0] = '\0' ;
    return (GrB_NOT_IMPLEMENTED) ;
    #endif
}

//------------------------------------------------------------------------------
// LAGraph_VectorTypeName: return the name of the GrB_Type of a GrB_Vector
//------------------------------------------------------------------------------

int LAGraph_VectorTypeName // returns 0 if successful, < 0 if failure
(
    char *name,             // name of the type of the vector v
    GrB_Vector v,           // vector to query
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (name != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // determine the name of the type of the GrB_Vector v
    //--------------------------------------------------------------------------

    #if LG_SUITESPARSE
    return (GxB_Vector_type_name (name, v)) ;
    #else
    name [0] = '\0' ;
    return (GrB_NOT_IMPLEMENTED) ;
    #endif
}

//------------------------------------------------------------------------------
// LAGraph_ScalarTypeName: return the name of the GrB_Type of a GrB_Scalar
//------------------------------------------------------------------------------

int LAGraph_ScalarTypeName // returns 0 if successful, < 0 if failure
(
    char *name,             // name of the type of the scalar v
    GrB_Scalar s,           // scalar to query
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (name != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // determine the name of the type of the GrB_Scalar s
    //--------------------------------------------------------------------------

    #if LG_SUITESPARSE
    return (GxB_Scalar_type_name (name, s)) ;
    #else
    name [0] = '\0' ;
    return (GrB_NOT_IMPLEMENTED) ;
    #endif
}

