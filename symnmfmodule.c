#include <Python.h>
#include "symnmf.h"
#include <assert.h>


/**
 * Converts a Python list of lists to a C array of doubles.
 *
 * @param points Python list of n sublists, each of length d
 * @param n Number of points (rows)
 * @param d Number of coordinates per point (columns)
 * @return Pointer to a newly allocated n×d array, or NULL on error
 */
double** extract_data(PyObject* points, int n, int d) {
    PyObject *sublist, *item;
    double **c_points;
    int i, j;

    // Allocate memory for the points array
    c_points =  create_mat(n, d);
    if (c_points == NULL) {
        return NULL;
    }

    // Populate the points array with values from the Python list
    for (i = 0; i < n; i++) {
        sublist = PyList_GetItem(points, i); // Get the i-th sublist
        if (!PyList_Check(sublist) || PyObject_Length(sublist) != d) {
            for (int k = 0; k <= i; k++) free(c_points[k]);
            free(c_points);
            return NULL;
        }

        for (j = 0; j < d; j++) {
            item = PyList_GetItem(sublist, j);
            if (PyFloat_Check(item)) {
                c_points[i][j] = PyFloat_AsDouble(item);
            } else if (PyLong_Check(item)) {
                c_points[i][j] = (double)PyLong_AsLong(item);
            } else {
                for (int k = 0; k <= i; k++) free(c_points[k]);
                free(c_points);
                return NULL;
            }
        }
    }
    return c_points;
}


/**
 * Converts a C 2D array to a Python list of lists.
 *
 * @param mat C array of size rows × cols
 * @param rows Number of rows
 * @param cols Number of columns
 * @return Python list of lists representing the matrix
 */
PyObject* matrix_to_pylist(double** mat, int rows, int cols) {
    PyObject *result;
    int i, j;
    result = PyList_New(rows);
    for (i = 0; i < rows; i++) {
        PyObject *centroid_row = PyList_New(cols);
        for (j = 0; j < cols; j++) {
            PyList_SetItem(centroid_row, j, PyFloat_FromDouble(mat[i][j]));
        }
        PyList_SetItem(result, i, centroid_row);
    }

    return result;
}


/**
 * Python wrapper for computing the symmetric matrix from input points.
 *
 * @param self Unused
 * @param args Tuple with one Python list of n points (each of length d)
 * @return Python list of lists representing the n×n symmetric matrix, or NULL on error
 */
static PyObject* SymWrapper(PyObject *self, PyObject *args)
{
    PyObject* points;      // Python list of lists representing points
    int n, d;              // number of points and dimension (inferred)
    double** c_points;     // C array of points
    double** A;            // the sym matrix to return
    PyObject* p_A;         // Python object to return

    if (!PyArg_ParseTuple(args, "O", &points)) {
        return NULL; // Error if input is invalid
    }

    n = PyList_Size(points);
    if (n == 0) return NULL;
    d = PyList_Size(PyList_GetItem(points, 0));

    // Use helper to extract points
    c_points = extract_data(points, n, d);
    if (c_points == NULL) {
        return NULL;
    }

    // Use helper to allocate A
    A = create_mat(n, n);
    if (A == NULL) {
        free_2D_array(c_points, n);
        return NULL;
    }

    // Fill the matrix
    get_sym(c_points, A, n, d);

    // Convert to Python list of tuples
    p_A = matrix_to_pylist(A, n, n);

    // Free allocated C memory
    free_2D_array(c_points, n);
    free_2D_array(A, n);
    return p_A;
}


/**
 * Python wrapper for computing the diagonal degree matrix.
 *
 * @param self Unused
 * @param args Tuple with one Python list of n points (each of length d)
 * @return Python list of lists representing the n×n diagonal degree matrix, or NULL on error
 */
static PyObject* DdgWrapper(PyObject *self, PyObject *args){
    PyObject* points;              // Python list of lists representing points
    int n, d;                      // number of points and dimension (inferred)
    double **c_points, **A, **D;       // C array of points, similarity matrix, diagonal matrix to return           
    PyObject* p_D;                 // Python object to return
    if (!PyArg_ParseTuple(args, "O", &points)) {
        return NULL;
    }
    n = PyList_Size(points);
    if (n == 0) return NULL;
    d = PyList_Size(PyList_GetItem(points, 0));
    
    c_points = extract_data(points, n, d);
    if (c_points == NULL) {
        return NULL;
    }

    A = create_mat(n, n);
    if (A == NULL){
        free_2D_array(c_points, n);
        return NULL;
    }

    D = create_mat(n, n);
    if (D == NULL){
        free_2D_array(c_points, n);
        free_2D_array(A, n);
        return NULL;
    }

    get_sym(c_points, A, n, d);
    get_diag(D, A, n);

    p_D = matrix_to_pylist(D, n, n);

    free_2D_array(c_points, n);
    free_2D_array(A, n);
    free_2D_array(D, n);
    return p_D;
}


/**
 * Python wrapper for computing the normalized matrix.
 *
 * @param self Unused
 * @param args Tuple with one Python list of n points (each of length d)
 * @return Python list of lists representing the n×n normalized matrix, or NULL on error
 */
static PyObject* NormWrapper(PyObject *self, PyObject *args){
    PyObject *points, *p_W; ;        // Python list of lists representing points, Python object to return
    int n, d;                       // number of points and dimension (inferred)
    double **c_points, **A, **D, **W;     // C array of points, similarity matrix, diagonal matrix, normalized matrix to return
    if (!PyArg_ParseTuple(args, "O", &points)) {
        return NULL;
    }
    n = PyList_Size(points);
    if (n == 0) return NULL;
    d = PyList_Size(PyList_GetItem(points, 0));
    c_points = extract_data(points, n, d);
    if (c_points == NULL) {
        return NULL;
    }
    A = create_mat(n, n);
    if (A == NULL){
        free_2D_array(c_points, n);
        return NULL;
    }
    D = create_mat(n, n);
    if (D == NULL){
        free_2D_array(c_points, n);
        free_2D_array(A, n);
        return NULL;
    }
    W = create_mat(n,n);
    if (W == NULL){
        free_2D_array(c_points, n);
        free_2D_array(A, n);
        free_2D_array(D, n);
        return NULL;
    }
    get_sym(c_points, A, n, d);
    get_diag(D, A, n);
    get_norm(W, D, A, n);
    p_W = matrix_to_pylist(W, n, n);
    free_2D_array(c_points, n);
    free_2D_array(A, n);
    free_2D_array(D, n);
    free_2D_array(W, n);
    return p_W;
}


/**
 * Python wrapper for computing the result of the symnmf algorithm.
 *
 * @param self Unused
 * @param args Tuple with two Python lists (W, H) and an integer k
 * @return Python list of lists representing the final H matrix, or NULL on error
 */
static PyObject* SymnmfWrapper(PyObject *self, PyObject *args){
    PyObject *py_W, *py_H, *p_final_H;
    int n, k, error_flag;
    double **W, **H, **final_H;    
    if (!PyArg_ParseTuple(args, "OOi", &py_W, &py_H, &k)) { 
        return NULL;
    }
    n = PyObject_Length(py_W);
    W = extract_data(py_W, n, n);                         
    if (W == NULL) {
        return NULL;
    }
    H = extract_data(py_H, n, k);
    if (H == NULL) {
        free_2D_array(W, n);
        return NULL;
    }
    final_H = create_mat(n, k);
    if (final_H == NULL){
        free_2D_array(W, n);
        free_2D_array(H, n);
        return NULL;
    }
    error_flag=get_symnmf(final_H, W, H, n, k);
    if (error_flag==1){
        free_2D_array(W, n);
        free_2D_array(H, n);
        free_2D_array(final_H, n);
        return NULL;
    }
    p_final_H = matrix_to_pylist(final_H, n, k);
    if (p_final_H == NULL) {                                    
        free_2D_array(W, n);
        free_2D_array(H, n);
        free_2D_array(final_H, n);
        return NULL;
    }
    free_2D_array(W, n);
    free_2D_array(H, n);
    free_2D_array(final_H, n);
    return p_final_H;
}


/**
 * Method table for Python bindings.
 * 
 * Includes functions for getting symmetric matrices, diagonal degree matrices, normalized matrices,
 * and running symnmf on W and H.
 */
static PyMethodDef PointsMethods[] = {
    {"sym", (PyCFunction)SymWrapper, METH_VARARGS, PyDoc_STR("get sym matrix")},
    {"ddg", (PyCFunction)DdgWrapper, METH_VARARGS, PyDoc_STR("get ddg matrix")},
    {"norm", (PyCFunction)NormWrapper, METH_VARARGS, PyDoc_STR("get norm matrix")},
    {"symnmf", (PyCFunction)SymnmfWrapper, METH_VARARGS, PyDoc_STR("run symnmf on W and H")},
    {NULL, NULL, 0, NULL} // Sentinel
};


/**
 * Module definition for the symnmfmodule.
 * 
 * Initializes the module with the symnmf-related functions.
 */
static struct PyModuleDef pointsmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule", // Module name
    NULL,            // Module documentation
    -1,              // Module state size
    PointsMethods    // Method table
};


/**
 * Initializes the symnmfmodule.
 * 
 * Creates and returns the module object for the symnmfmodule.
 * Returns NULL on failure.
 */
PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    PyObject* mod;
    mod = PyModule_Create(&pointsmodule);
    if (!mod)
    {
        return NULL;
    }
    return mod;
}
