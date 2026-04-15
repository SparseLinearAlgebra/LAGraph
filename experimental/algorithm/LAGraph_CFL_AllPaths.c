#define LG_FREE_WORK          \
  {                           \
    GrB_free(&AllPaths_semiring);   \
    GrB_free(&AllPaths_monoid);     \
    GrB_free(&AllPaths_monoid_get_nvals);     \
    GrB_free(&bottom_scalar); \
    GrB_free(&IAllPaths_mult);      \
    GrB_free(&AllPaths_set);        \
    GrB_free(&AllPaths_mult);       \
    GrB_free(&AllPaths_add);        \
    GrB_free(&AllPaths_add_get_nvals);        \
    GrB_free(&Theta);         \
  }

#include "LG_internal.h"
#include <LAGraphX.h>

//Merging two ordered arrays of internal vertices in the add function
static GrB_Index* merge_all_paths(GrB_Index* n, const GrB_Index* a, const size_t na, const GrB_Index* b, const size_t nb)
{  
  GrB_Index *tmp = malloc((na + nb) * sizeof(GrB_Index));
  
  size_t ia = 0, ib = 0, outn = 0;
  // Handle the first elements of both arrays to initialize tmp and avoid checking if tmp is empty in the loop
  if (na > 0 && nb > 0) {
    if (a[0] < b[0]) {
      tmp[outn++] = a[ia++];
    } else if (b[0] < a[0]) {
      tmp[outn++] = b[ib++];
    } else {
      tmp[outn++] = a[ia++];
      ib++;
    }
  } else if (na > 0) {
    tmp[outn++] = a[ia++];
  } else {
    tmp[outn++] = b[ib++];
  }
  
  while (ia < na && ib < nb) {
    GrB_Index va = a[ia];
    GrB_Index vb = b[ib];
    if (va < vb) {
      if (tmp[outn-1] != va) tmp[outn++] = va;
      ia++;
    } else if (vb < va) {
      if (tmp[outn-1] != vb) tmp[outn++] = vb;
      ib++;
    } else {
      if (tmp[outn-1] != va) tmp[outn++] = va;
      ia++; ib++;
    }
  }
  while (ia < na) {
    GrB_Index va = a[ia++];
    if (tmp[outn-1] != va) tmp[outn++] = va;
  }
  while (ib < nb) {
    GrB_Index vb = b[ib++];
    if (tmp[outn-1] != vb) tmp[outn++] = vb;
  }
  
  tmp = realloc(tmp, outn * sizeof(GrB_Index));
  *n = outn;
  
  return tmp;
}

static GrB_Index* insert_all_paths(size_t* n, const GrB_Index* arr, size_t len, GrB_Index value)
{
  // Array is ordered, so we can use binary search to find the position to insert value
  size_t l = 0, r = len, m = 0;
  while (l < r) {
    m = l + (r - l) / 2;
    if (arr[m] < value) l = m + 1;
    else r = m;
  }

    // If value is already in the array, return the original array
    if (l < len && arr[l] == value) {
        GrB_Index *tmp = malloc(len * sizeof(GrB_Index));
        memcpy(tmp, arr, len * sizeof(GrB_Index));
        *n = len;
        return tmp;
    }

    GrB_Index *tmp = malloc((len + 1) * sizeof(GrB_Index));
    memcpy(tmp, arr, l * sizeof(GrB_Index));
    tmp[l] = value;
    memcpy(tmp + l + 1, arr + l, (len - l) * sizeof(GrB_Index));

    *n = len + 1;
    return tmp;
}

void add_all_paths(AllPathsElem *z, AllPathsElem *x, AllPathsElem *y)
{
    // temp is needed to avoid freeing the memory of z in case z == x or z == y
    AllPathsElem temp;

    // x->n and y->n are not equal to zero
    if (x->n == 1 && y->n == 1) {
        if (x->data.single_elem == y->data.single_elem) {
            temp.n = 1;
            temp.data.single_elem = x->data.single_elem;
        } else {
            temp.n = 2;
            temp.data.middle = malloc(2 * sizeof(GrB_Index));
            if (x->data.single_elem < y->data.single_elem) {
                temp.data.middle[0] = x->data.single_elem;
                temp.data.middle[1] = y->data.single_elem;
            } else {
                temp.data.middle[0] = y->data.single_elem;
                temp.data.middle[1] = x->data.single_elem;
            }
        }
    } else if (x->n == 1) {
        temp.data.middle = insert_all_paths(&temp.n, y->data.middle, y->n, x->data.single_elem);
    } else if (y->n == 1) {
        temp.data.middle = insert_all_paths(&temp.n, x->data.middle, x->n, y->data.single_elem);
    } else {
        temp.data.middle = merge_all_paths(&temp.n, x->data.middle, x->n, y->data.middle, y->n);
    }

    if(x->n > 1){
      free(x->data.middle);
    }
    if(y->n > 1){
      free(y->data.middle);
    }

    *z = temp;
}

void mult_all_paths(AllPathsElem *z,
                     const AllPathsElem *x, GrB_Index ix, GrB_Index jx,
                     const AllPathsElem *y, GrB_Index iy, GrB_Index jy,
                     const void *theta)
{
  z->data.single_elem = jx;
  z->n = 1;
}

void set_all_paths(AllPathsElem *z, const AllPathsElem *x, const bool *edge_exist)
{
  z->data.single_elem = GrB_INDEX_MAX; // A special value to indicate that this path corresponds to a terminal rule (A->t) or an epsilon rule (A->eps)
  z->n = 1;
}

#define MULT_PATH_INDEX_DEFN                                                   \
"void mult_all_paths(AllPathsElem *z, \n"      \
"                     const AllPathsElem *x, GrB_Index ix, GrB_Index jx, \n"      \
"                     const AllPathsElem *y, GrB_Index iy, GrB_Index jy, \n"      \
"                     const void *theta) \n"      \
"{ \n"      \
"  z->data.single_elem = jx; \n"      \
"  z->n = 1; \n"      \
"}"

//Adding the count of all internal vertices in a reduction
void add_get_nvals_all_paths(AllPathsElem *z, const AllPathsElem *x, const AllPathsElem *y)
{
  z->n = x->n + y->n;
}

//Made global so that the get_nvals_all_paths matches the GrB_Matrix_nvals signature
static GrB_Type AllPaths_type = NULL;
static GrB_Monoid AllPaths_monoid_get_nvals = NULL;

//A function that replaces GrB_Matrix_nvals in Reachability to check if new vertices have been added to the matrix.
static GrB_Info get_nvals_all_paths(GrB_Index *nvals, const GrB_Matrix A)
{
  GrB_Scalar s = NULL;
  GrB_Scalar_new(&s, AllPaths_type);
  GrB_Info info = GrB_reduce(s, NULL, AllPaths_monoid_get_nvals, A, NULL);
  
  if (info != GrB_SUCCESS)
  {
    GrB_free(&s);
    return info;
  }
  
  AllPathsElem result = {0};
  info = GrB_Scalar_extractElement_UDT(&result, s);
  if (info == GrB_NO_VALUE)
  {
      result.n = 0;
  }
  else if (info != GrB_SUCCESS)
  {
      GrB_free(&s);
      return info;
  }
    
  *nvals = result.n;
  GrB_free(&s);
  return GrB_SUCCESS;
}

GrB_Info LAGraph_CFL_AllPaths(
    // Output
    GrB_Matrix *outputs, // Array of matrices containing results.
                         // The size of the array must be equal to nonterms_count.
                         // Matrix elements are ordered arrays of intermediate vertices type AllPathsElem.
                         // Before free outputs[k], you need to free arrays from all matrix elements.
                         // For all values of M from the array of the matrix element
                         // outputs[k] on the I row of the J column:
                         // There are paths from I to M by nonterminal N1 and from M to J by nonterminal N2,
                         // and A->N1 N2 where outputs[k] corresponds to nonterminal A.
                         // GrB_INDEX_MAX in the array is a special value for A->eps and A->t.
    // AllPaths type - elements of the output matrices.
    GrB_Type *all_paths_ptr_t,      // Pass a pointer to GrB_Type.
    // Input
    const GrB_Matrix *adj_matrices, // Array of adjacency matrices representing the graph.
                                    // The length of this array is equal to the count of
                                    // terminals (terms_count).
                                    //
                                    // adj_matrices[t]: (i, j) == 1 if and only if there
                                    // is an edge between nodes i and j with the label of
                                    // the terminal corresponding to index 't' (where t is
                                    // in the range [0, terms_count - 1]).
    int64_t terms_count,            // The total number of terminal symbols in the CFG.
    int64_t nonterms_count,         // The total number of non-terminal symbols in the CFG.
    const LAGraph_rule_WCNF *rules, // The rules of the CFG.
    int64_t rules_count,            // The total number of rules in the CFG.
    char *msg                       // Message string for error reporting.
)
{
#if GxB_IMPLEMENTATION < GxB_VERSION(9, 4, 5)
  return (GrB_NOT_IMPLEMENTED);
#else
  // Semiring components
  GrB_BinaryOp AllPaths_add = NULL;
  GrB_BinaryOp AllPaths_add_get_nvals = NULL;
  GrB_Monoid AllPaths_monoid = NULL;
  GxB_IndexBinaryOp IAllPaths_mult = NULL;
  GrB_BinaryOp AllPaths_mult = NULL;
  GrB_Semiring AllPaths_semiring = NULL;
  GrB_BinaryOp AllPaths_set = NULL;
  GrB_Scalar Theta = NULL;
  GrB_Scalar bottom_scalar = NULL;
  
  GrB_free(all_paths_ptr_t);
  GRB_TRY(GrB_Type_new(all_paths_ptr_t, sizeof(AllPathsElem)));
  AllPaths_type = *all_paths_ptr_t;
  
  GRB_TRY(GrB_Scalar_new(&Theta, GrB_BOOL));
  GRB_TRY(GrB_Scalar_setElement_BOOL(Theta, false));
  
  AllPathsElem bottom = {0};
  GRB_TRY(GrB_Scalar_new(&bottom_scalar, AllPaths_type));
  GRB_TRY(GrB_Scalar_setElement_UDT(bottom_scalar, (void *)(&bottom)));
  
  GRB_TRY(GrB_BinaryOp_new(
                           &AllPaths_add,
                           (void *)add_all_paths,
                           AllPaths_type,
                           AllPaths_type,
                           AllPaths_type));
  
  GRB_TRY(GrB_Monoid_new(
                         &AllPaths_monoid,
                         AllPaths_add,
                         (void *)(&bottom)));
  
  GRB_TRY(GxB_IndexBinaryOp_new(
                                &IAllPaths_mult,
                                (void *)mult_all_paths,
                                AllPaths_type,
                                AllPaths_type,
                                AllPaths_type,
                                GrB_BOOL,
                                "mult_all_paths",
                                MULT_PATH_INDEX_DEFN));
  
  GRB_TRY(GxB_BinaryOp_new_IndexOp(
                                   &AllPaths_mult,
                                   IAllPaths_mult,
                                   Theta));
  
  GRB_TRY(GrB_Semiring_new(
                           &AllPaths_semiring,
                           AllPaths_monoid,
                           AllPaths_mult));
  
  GRB_TRY(GrB_BinaryOp_new(
                           &AllPaths_set,
                           (void *)set_all_paths,
                           AllPaths_type,
                           AllPaths_type,
                           GrB_BOOL));
  
  GRB_TRY(GrB_BinaryOp_new(
                           &AllPaths_add_get_nvals,
                           (void *)add_get_nvals_all_paths,
                           AllPaths_type,
                           AllPaths_type,
                           AllPaths_type));
  
  GRB_TRY(GrB_Monoid_new(
                         &AllPaths_monoid_get_nvals,
                         AllPaths_add_get_nvals,
                         (void *)(&bottom)));
  
  CFL_Semiring semiring = {
      .type = AllPaths_type,
      .semiring = AllPaths_semiring,
      .add = AllPaths_add,
      .mult = AllPaths_mult,
      .init_path = AllPaths_set,
      .bottom_scalar = bottom_scalar,
      .get_nvals = get_nvals_all_paths};
  
  LG_TRY(LAGraph_CFPQ_core(outputs, adj_matrices, terms_count, nonterms_count, rules, rules_count, &semiring, msg));
  
  LG_FREE_WORK;
  return GrB_SUCCESS;
#endif
}

// Helper function to free the output matrix of LAGraph_CFL_AllPaths, which contains elements of type AllPathsElem with dynamically allocated arrays of intermediate vertices. 
static void free_AllPaths_matrix(GrB_Matrix* ptr_output) 
{
  GxB_Iterator iterator;
  GxB_Iterator_new(&iterator);
  GrB_Info info = GxB_Matrix_Iterator_attach(iterator, *ptr_output, NULL);
  info = GxB_Matrix_Iterator_seek(iterator, 0);
  AllPathsElem val;

  while (info != GxB_EXHAUSTED)
  {
    GxB_Iterator_get_UDT(iterator, (void*)&val);
    if (val.n > 1 && val.data.middle != NULL) {
      free(val.data.middle);
    }
    info = GxB_Matrix_Iterator_next(iterator);
  }

  GrB_free(&iterator);
  GrB_free(ptr_output);
}

// Free outputs and all_paths_ptr_t after you have finished working with the output matrices from LAGraph_CFL_AllPaths.
// do outputs = NULL, all_paths_ptr_t = NULL after LAGraph_CFL_AllPaths_free_outputs
GrB_Info LAGraph_CFL_AllPaths_free_outputs(GrB_Matrix* outputs, int64_t nonterms_count, GrB_Type* all_paths_ptr_t)
{
#if GxB_IMPLEMENTATION < GxB_VERSION(9, 4, 5)
  return (GrB_NOT_IMPLEMENTED);
#else
  if (outputs) {
    for (size_t i = 0; i < nonterms_count; i++) {
      if (outputs[i] == NULL)
        continue;
      free_AllPaths_matrix(&outputs[i]);
      outputs[i] = NULL;
    }
    free(outputs);
  }
  GrB_free(all_paths_ptr_t);
  return GrB_SUCCESS;
#endif
}

