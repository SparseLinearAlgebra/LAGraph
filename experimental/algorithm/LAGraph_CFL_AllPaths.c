#define LG_FREE_WORK          \
  {                           \
    GrB_free(&AllPaths_semiring_free);   \
    GrB_free(&AllPaths_monoid_free);     \
    GrB_free(&bottom_scalar); \
    GrB_free(&IAllPaths_mult);      \
    GrB_free(&AllPaths_set);        \
    GrB_free(&AllPaths_mult);       \
    GrB_free(&AllPaths_add);        \
    GrB_free(&AllPaths_add_free);        \
    GrB_free(&Theta);         \
  }

#include "LG_internal.h"
#include <LAGraphX.h>

static GrB_Index* merge_all_paths(GrB_Index* n, const void* left, const GrB_Index na, const void* right, const GrB_Index nb){
  GrB_Index* a = (GrB_Index*) left;
  GrB_Index* b = (GrB_Index*) right;
  
  if (na == 0 && nb == 0) {
      *n = 0;
      return NULL;
  }
  
  GrB_Index *tmp = malloc((na + nb) * sizeof(GrB_Index));
  //    LG_TRY(LAGraph_Malloc((void**)&tmp, na+nb, sizeof(GrB_Index), msg));
  
  GrB_Index ia = 0, ib = 0, outn = 0;
  while (ia < na && ib < nb) {
    GrB_Index va = a[ia];
    GrB_Index vb = b[ib];
    if (va < vb) {
      if (outn == 0 || tmp[outn-1] != va) tmp[outn++] = va;
      ia++;
    } else if (vb < va) {
      if (outn == 0 || tmp[outn-1] != vb) tmp[outn++] = vb;
      ib++;
    } else {
      if (outn == 0 || tmp[outn-1] != va) tmp[outn++] = va;
      ia++; ib++;
    }
  }
  while (ia < na) {
    GrB_Index va = a[ia++];
    if (outn == 0 || tmp[outn-1] != va) tmp[outn++] = va;
  }
  while (ib < nb) {
    GrB_Index vb = b[ib++];
    if (outn == 0 || tmp[outn-1] != vb) tmp[outn++] = vb;
  }
  
  //    LG_TRY(LAGraph_Realloc((void**)&tmp, outn, na+nb, sizeof(GrB_Index), msg));
  GrB_Index *sh = realloc(tmp, outn * sizeof(GrB_Index));
  if (sh) tmp = sh;
  
  *n = outn;
  return tmp;
}

void clear_all_paths_vex(AllPathsVex *z){
//  if(z->middle) LG_TRY(LAGraph_Free((void**) z->middle, msg));
  if(z->middle){
    free(z->middle);
    z->middle=NULL;
  }
  z->n = 0;
}

void add_all_paths(AllPathsVex *z, const AllPathsVex *x, const AllPathsVex *y)
{
  z->middle = merge_all_paths(&z->n, x->middle, x->n, y->middle, y->n);
}

void add_all_paths_free(AllPathsVex *z, AllPathsVex *x, AllPathsVex *y)
{
  AllPathsVex v_temp;
  AllPathsVex* temp = &v_temp;
  temp->middle = merge_all_paths(&temp->n, x->middle, x->n, y->middle, y->n);
  clear_all_paths_vex(x);
  clear_all_paths_vex(y);
//  clear_all_paths_vex(z);
  z->middle = temp->middle;
  z->n = temp->n;
}

void mult_all_paths_index(AllPathsVex *z,
                     const AllPathsVex *x, GrB_Index ix, GrB_Index jx,
                     const AllPathsVex *y, GrB_Index iy, GrB_Index jy,
                     const void *theta)
{
//  LG_TRY(LAGraph_Malloc((void**) &z->middle, 1, size_of(GrB_Index), msg));
  z->middle = malloc(sizeof(GrB_Index));
  z->middle[0] = jx;
  z->n = 1;
}

void set_all_paths(AllPathsVex *z, const AllPathsVex *x, const bool *edge_exist)
{
  AllPathsVex temp;
  temp.middle = NULL;
  temp.n = 0;
  if (edge_exist && *edge_exist){
    temp.middle = malloc(sizeof(GrB_Index));
    temp.n = 1;
    temp.middle[0] = GrB_INDEX_MAX;
  }
  else{
    temp.middle = NULL;
    temp.n = 0;
  }
  z->middle = merge_all_paths(&z->n, x->middle, x->n, temp.middle, temp.n);
  clear_all_paths_vex(&temp);
}

#define MULT_PATH_INDEX_DEFN                                                   \
"void mult_all_paths_index(AllPathsVex *z, \n"      \
"                     const AllPathsVex *x, GrB_Index ix, GrB_Index jx, \n"      \
"                     const AllPathsVex *y, GrB_Index iy, GrB_Index jy, \n"      \
"                     const void *theta) \n"      \
"{ \n"      \
"  //clear_all_paths_vex(z); \n"      \
"  z->middle = malloc(sizeof(GrB_Index)); \n"      \
"  z->middle[0] = jx; \n"      \
"  z->n = 1; \n"      \
"}"

//A function that replaces GrB_Matrix_nvals for counting non-zero elements
GrB_Info all_paths_get_nvals(GrB_Index *nvals, const GrB_Matrix A){
  GrB_Index nnz = 0;
  GrB_Index accum = 0;
  GrB_Matrix_nvals(&nnz, A);
  *nvals = 0;
  void *val_void = NULL;
  
  if(nnz==0){
    *nvals = 0;
    return GrB_SUCCESS;
  }
  
  val_void = malloc(nnz * sizeof(AllPathsVex));

  GrB_Matrix_extractTuples(NULL, NULL, val_void, &nnz, A);
  
  AllPathsVex *val = (AllPathsVex *) val_void;
  
  size_t end = nnz;
  
  for (size_t start = 0; start < end; ++start){
    accum+=val[start].n;
  }
  
  *nvals = accum;
  
  free(val_void);
  
  return GrB_SUCCESS;
}

GrB_Info LAGraph_CFL_AllPaths(
    // Output
    GrB_Matrix *outputs, // Array of matrices containing results.
                         // The size of the array must be equal to nonterms_count.
                         //
                         // outputs[k]: (i, j) contains a AllPathsVex structure if and only if there is a path
                         // from node i to node j whose edge labels form a word
                         // derivable from the non-terminal 'k' of the specified CFG.
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
  // Semiring components
  GrB_Type AllPaths_type = NULL;
  GrB_BinaryOp AllPaths_add = NULL;
  GrB_BinaryOp AllPaths_add_free = NULL;
  GrB_Monoid AllPaths_monoid_free = NULL;
  GxB_IndexBinaryOp IAllPaths_mult = NULL;
  GrB_BinaryOp AllPaths_mult = NULL;
  GrB_Semiring AllPaths_semiring_free = NULL;
//  GxB_IndexBinaryOp IAllPaths_set = NULL;
  GrB_BinaryOp AllPaths_set = NULL;
  GrB_Scalar Theta = NULL;
  GrB_Scalar bottom_scalar = NULL;

  GRB_TRY(GrB_Type_new(&AllPaths_type, sizeof(AllPathsVex)));

  GRB_TRY(GrB_Scalar_new(&Theta, GrB_BOOL));
  GRB_TRY(GrB_Scalar_setElement_BOOL(Theta, false));

  AllPathsVex bottom = {0, NULL};
  GRB_TRY(GrB_Scalar_new(&bottom_scalar, AllPaths_type));
  GRB_TRY(GrB_Scalar_setElement_UDT(bottom_scalar, (void *)(&bottom)));

  GRB_TRY(GrB_BinaryOp_new(
      &AllPaths_add,
      (void *)add_all_paths,
      AllPaths_type,
      AllPaths_type,
      AllPaths_type));
  
  GRB_TRY(GrB_BinaryOp_new(
      &AllPaths_add_free,
      (void *)add_all_paths_free,
      AllPaths_type,
      AllPaths_type,
      AllPaths_type));

  GRB_TRY(GrB_Monoid_new(
      &AllPaths_monoid_free,
      AllPaths_add_free,
      (void *)(&bottom)));
  
  GRB_TRY(GxB_IndexBinaryOp_new(
      &IAllPaths_mult,
      (void *)mult_all_paths_index,
      AllPaths_type,
      AllPaths_type,
      AllPaths_type,
      GrB_BOOL,
      "mult_all_paths_index",
      MULT_PATH_INDEX_DEFN));

  GRB_TRY(GxB_BinaryOp_new_IndexOp(
      &AllPaths_mult,
      IAllPaths_mult,
      Theta));

  GRB_TRY(GrB_Semiring_new(
      &AllPaths_semiring_free,
      AllPaths_monoid_free,
      AllPaths_mult));
  
  GRB_TRY(GrB_BinaryOp_new(
                           &AllPaths_set,
                           (void *)set_all_paths,
                            AllPaths_type,
                            AllPaths_type,
                            GrB_BOOL));
  
//  GRB_TRY(GxB_IndexBinaryOp_new(
//      &IAllPaths_set,
//      (void *)set_all_paths_index,
//      AllPaths_type,
//      AllPaths_type,
//      GrB_BOOL,
//      GrB_BOOL,
//      "set_all_paths_index",
//      SET_PATH_INDEX_DEFN));

//  GRB_TRY(GxB_BinaryOp_new_IndexOp(
//      &AllPaths_set,
//      IAllPaths_set,
//      Theta));
  
  CFL_Semiring semiring = {.type = AllPaths_type,
      .semiring = AllPaths_semiring_free,
      .add = AllPaths_add_free,
      .add_eps = AllPaths_add,
      .mult = AllPaths_mult,
      .init_path = AllPaths_set,
      .bottom_scalar = bottom_scalar,
      .get_nvals = all_paths_get_nvals};
  
  LG_TRY(LAGraph_CFPQ_core(outputs, adj_matrices, terms_count, nonterms_count, rules, rules_count, &semiring, msg));
  LG_FREE_WORK;
  return GrB_SUCCESS;
}
