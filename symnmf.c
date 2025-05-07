#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "symnmf.h"

#define MAX_ITER 300
#define EPSILON 1e-4

const char *CASE_ERROR =  "An Error Has Occurred\n";

int run_program(double** points, char* goal, int N, int d);
int check_input(int argc, char *argv[]);
int get_dimensions(char *path, int *N, int *D);
int points_file_to_mat(char *path, double **data, int *N, int *d);
void free_2D_array(double **array, int rows);
double get_distance(double *p1, double *p2, int d);
double** create_mat(int num_of_rows, int num_of_cols);
void print_mat(double **mat, int rows, int cols);
void get_sym(double **points, double **A, int N, int d);
void get_diag(double **D, double **A, int N);
void get_norm(double **W, double **D, double **A, int N);
void matrix_multiply(double **C,double **A, double **B, int m, int k, int n);
void matrix_transpose(double **T,double **A, int rows, int cols);
double Frobenius_norm(double** A, double** B, int rows, int cols);
void copy_matrix(double** source, double** destination, int n, int k);
int make_step(double** new_H ,double** W,double** H_just_before_the_step, int n, int k);
int get_symnmf(double** new_H, double **W, double **H, int n, int k);

/**
 * Main function to process input arguments, read the dataset, and run the specified matrix computation.
 *
 * The function:
 * - Validates the input arguments.
 * - Reads the dataset from a file and loads it into a matrix.
 * - Performs matrix computations based on the specified goal (e.g., "sym", "ddg", or "norm").
 * - Frees any allocated memory and handles errors during execution.
 *
 * Parameters:
 *   argc - Number of command-line arguments.
 *   argv - Array of command-line arguments.
 *
 * Returns:
 *   0 if the program runs successfully, 1 if any error occurs (e.g., invalid input, file reading failure).
 */
int main(int argc, char *argv[]) {
    double **points;
    int d, N, error_flag;
    char *goal = NULL, *path = NULL;

    error_flag = check_input(argc, argv);
    if (error_flag==1) {
        printf("%s", CASE_ERROR);
        return 1;
    }
    goal = argv[1];
    path = argv[2];
    
    error_flag=get_dimensions(path, &N, &d);
    if (error_flag==1) {
        printf("%s", CASE_ERROR);
        return 1;
    }

    points = create_mat(N, d); /*allocate matrix*/
    if (points == NULL) {
        printf("%s", CASE_ERROR);
        return 1;
    }
 
    error_flag = points_file_to_mat(path, points, &N, &d); /*get data from file to the matrix*/
    if (error_flag==1) {
        free_2D_array(points, N);
        printf("%s", CASE_ERROR);
        return 1;
    }

    error_flag = run_program(points, goal, N, d);
    free_2D_array(points, N);

    if (error_flag==1) {
        printf("%s", CASE_ERROR);
        return 1;
    }
    return 0;
}


/**
 * Validates the command-line arguments for the program.
 *
 * The function checks:
 * - That the correct number of arguments is provided (exactly 3).
 * - That the goal argument is one of: "sym", "ddg", or "norm".
 * - That the file path ends with ".txt".
 *
 * Parameters:
 *   argc - Number of command-line arguments.
 *   argv - Array of argument strings.
 *
 * Returns:
 *   0 if all inputs are valid, 1 otherwise.
 */
int check_input(int argc, char *argv[])
{
    char *goal, *path;
    size_t len;

    if (argc != 3) {
        return 1;
    }
    goal = argv[1];
    path = argv[2];

    if ((strcmp(goal, "sym") != 0) && (strcmp(goal, "ddg") != 0) && (strcmp(goal, "norm") != 0)) {
        return 1;
    }

    len = strlen(path);

    if (len < 4 || strcmp(path + len - 4, ".txt") != 0) {
        return 1;
    }
    return 0;
}


/**
 * Runs the matrix computation based on the specified goal.
 *
 * This function:
 * - Creates a matrix `A` and computes a symmetric matrix based on `points`.
 * - If the goal is "sym", it prints the matrix `A`.
 * - If the goal is "ddg" or "norm", it computes a diagonal matrix `D` and:
 *   - For "ddg", prints `D`.
 *   - For "norm", computes a normalized matrix `W`, prints it, and frees allocated memory.
 *
 * Parameters:
 *   points - A 2D array representing the input data points.
 *   goal - A string indicating the goal of the computation ("sym", "ddg", or "norm").
 *   N - The number of rows/columns in the matrices.
 *   d - The dimensionality of the data points.
 *
 * Returns:
 *   0 if the program ran successfully, 1 if there was an error (e.g., memory allocation failure).
 */
int run_program(double** points, char* goal, int N, int d){
    double **A, **W, **D;
    A = create_mat(N,N);
    if (A == NULL) {
        return 1;
    }
    get_sym(points, A, N, d);
        if (strcmp(goal, "sym") == 0) {
            print_mat(A,N,N);
        }

    if (strcmp(goal, "ddg") == 0 || strcmp(goal, "norm") == 0) {
        D = create_mat(N,N);

    if (D == NULL) {
        free_2D_array(A, N);
        return 1;
    }
        get_diag(D,A,N);
        if (strcmp(goal, "ddg") == 0) {
            print_mat(D,N,N);
            free_2D_array(D, N);
        }
    } 
    if (strcmp(goal, "norm") == 0){
        W = create_mat(N,N);

        if (W == NULL) {
            free_2D_array(D, N);
            free_2D_array(A, N);
            return 1;
        }
        get_norm(W,D,A,N);
        print_mat(W,N,N);
        free_2D_array(D, N);
        free_2D_array(W, N);
    }
    free_2D_array(A, N);
    return 0;
}


/*
 * FUNCTION: get_dimensions
 * -------------------------
 * Reads a file containing a dataset and assigns the number of data points (N) 
 * and the number of dimensions (D) of each point.
 *
 * INPUTS:
 *   - char *path: Path to the input file (CSV format expected).
 *   - int *N: Pointer to an integer to store the number of rows (data points).
 *   - int *D: Pointer to an integer to store the number of dimensions per point.
 *
 * OUTPUT:
 *   - Returns 0 if successful, or 1 if an error occurs (e.g., file not found).
 *   - Updates *N with the number of data points (rows).
 *   - Updates *D with the number of dimensions per data point (columns).
 */
 
int get_dimensions(char *path, int *N, int *D) {
    int c, prev = 0;

    /* Open the file for reading */
    FILE *file = fopen(path, "r");
    if (file == NULL) {
        return 1; /* File not found, return error */
    }

    /* Initialize dimensions and rows to zero */
    *D = 0;
    *N = 0;

    /* Count dimensions in the first line */
    while ((c = fgetc(file)) != '\n' && c != EOF) {
        if (c == ',') {
            (*D)++;  /* Each comma indicates an additional dimension */
        }
        prev = c;
    }
    /* If there are no commas in the first line, we assume there's one dimension */
    if (*D == 0 && prev != '\n') {
        *D = 1;  /* If no commas, assume single column (k=1) */
    }
    (*D)++;  /* Account for the last value after the last comma or for the only value in case of k=1 */

    /*Count the rows in the file */
    fseek(file, 0, SEEK_SET);  /* Reset file pointer to the beginning */
    while ((c = fgetc(file)) != EOF) {
        if (c == '\n' || c == '\r') {  /* Handle both '\n' and '\r' line endings */
            (*N)++;  /* Increment row count at each newline or carriage return */
        }
        prev = c;
    }
    /* If the last line doesn't end with a newline, increment row count */
    if (prev != '\n' && prev != '\r' && prev != EOF) {
        (*N)++;
    }
    fclose(file);  /* Close the file after counting rows and dimensions */
    return 0;  /* Return success */
}



/*
 * FUNCTION: points_file_to_mat
 * ----------------------------
 * Reads a CSV file and populates a pre-allocated matrix with data points.
 *
 * INPUTS:
 *   - char *path: Path to the input file (CSV format expected).
 *   - double **data: Pointer to a pre-allocated 2D array (N x d) to store the data points.
 *   - int *N: Pointer to an integer indicating the expected number of data points (rows).
 *   - int *d: Pointer to an integer indicating the number of dimensions per point (columns).
 *
 * OUTPUT:
 *   - Returns 0 if successful, or 1 if an error occurs (e.g., file not found).
 *   - Fills the provided data matrix with the parsed values from the file.
 */
int points_file_to_mat(char *path, double **data, int *N, int *d){
    FILE *file;
    char line[1024];
    int i = 0, j = 0;
    /* Read data */
    file = fopen(path, "r");
    if (file == NULL) {
        return 1;
    }

    while (fgets(line, sizeof(line), file) != NULL && i < *N) {
        char *ptr;
        data[i][0] = strtod(line, &ptr);
        for (j=1; j < *d; j++) {
            data[i][j] = strtod(ptr + 1, &ptr);
        }
        i++;
    }
    fclose(file);
    return 0;
}


/*
 * FUNCTION: free_2D_array
 * -----------------------
 * Frees memory allocated for a 2D array.
 *
 * INPUTS:
 *   - double **array: Pointer to the array of pointers (2D array) to be freed.
 *   - int rows: Number of rows (i.e., number of pointers in the outer array).
 *
 * OUTPUT:
 *   - None. Frees all allocated memory for the 2D array.
 */
void free_2D_array(double **array, int rows) {
    int i;
    for (i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
}


/*
 * FUNCTION: get_distance
 * ----------------------
 * Computes the squared Euclidean distance between two points in d-dimensional space.
 *
 * INPUTS:
 *   - double *p1: Pointer to the first point (array of doubles).
 *   - double *p2: Pointer to the second point (array of doubles).
 *   - int d: Number of dimensions of the points.
 *
 * OUTPUT:
 *   - Returns the squared Euclidean distance between the two points.
 */
double get_distance(double *p1, double *p2, int d) {
    int i;
    double delta = 0,dist = 0;
    for(i=0;i<d;i++){
        delta = p1[i]-p2[i];
        dist += delta*delta;
    }
    return dist;
}


/*
 * FUNCTION: create_mat
 * --------------------
 * Allocates and initializes a 2D matrix of doubles with all elements set to zero.
 *
 * INPUTS:
 *   - int num_of_rows: Number of rows in the matrix.
 *   - int num_of_cols: Number of columns in the matrix.
 *
 * OUTPUT:
 *   - Returns a pointer to the newly allocated 2D matrix, or NULL if allocation fails.
 */
double **create_mat(int num_of_rows, int num_of_cols) {
    double **mat;
    int i,j;
    mat = (double**)calloc(num_of_rows, sizeof(double*));
    if (mat==NULL) {
        return NULL;
    }
    for (i = 0; i < num_of_rows; i++) {
        mat[i] =  (double*)calloc(num_of_cols, sizeof(double));
        if (mat[i] == NULL) { 
            for (j = 0; j < i; j++) {
                free(mat[j]);
            }
            free(mat);
            return NULL; 
        }
    }
    return mat;
}


/*
 * FUNCTION: print_mat
 * -------------------
 * Prints the contents of a 2D matrix of doubles with 4 decimal places of precision.
 *
 * INPUTS:
 *   - double **mat: Pointer to the 2D matrix to print.
 *   - int rows: Number of rows in the matrix.
 *   - int cols: Number of columns in the matrix.
 *
 * OUTPUT:
 *   - None. Prints the matrix to standard output.
 */
void print_mat(double **mat, int rows, int cols) {
    int i, j;
    
    for (i = 0; i < rows; i++) {  
        for (j = 0; j < cols; j++) {  
            if (j > 0) {
                printf(",");  
                }
            printf("%.4f", mat[i][j]); 
            }
        printf("\n");  
        }
    }


/*
 * FUNCTION: get_sym
 * -----------------
 * Computes a symmetric matrix A based on pairwise distances between points.
 *
 * INPUTS:
 *   - double **points: Pointer to the dataset (N points, each of dimension d).
 *   - double **A: Pointer to the matrix to store the result (N x N), modified in place.
 *   - int N: Number of points.
 *   - int d: Dimension of each point.
 *
 * OUTPUT:
 *   - None. The matrix A is updated in place with computed values.
 */
void get_sym(double **points, double **A, int N, int d) {
    int i, j;
    double dist;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i < j) {
                dist = get_distance(points[i], points[j], d);
                A[i][j] = exp(-dist / 2);
            }
            if (i > j) {
                A[i][j] = A[j][i];
            }
        }
    }
}


/*
 * FUNCTION: get_diag
 * ------------------
 * Computes the diagonal degree matrix D from a symmetric matrix A.
 *
 * INPUTS:
 *   - double **D: Pointer to the matrix to store the result (N x N), modified in place.
 *   - double **A: Pointer to the symmetric matrix (N x N).
 *   - int N: Dimension of the matrices.
 *
 * OUTPUT:
 *   - None. The matrix D is updated in place as a diagonal degree matrix.
 */
void get_diag(double **D, double **A, int N) {
    int i,j;
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            D[i][i]+=A[i][j];
        }
    }
}


/*
 * FUNCTION: get_norm
 * ------------------
 * Computes the normalized matrix W using the symmetric matrix A and diagonal matrix D.
 *
 * INPUTS:
 *   - double **W: Pointer to the matrix to store the result (N x N), modified in place.
 *   - double **D: Pointer to the diagonal degree matrix (N x N), modified in place to D^(-1/2).
 *   - double **A: Pointer to the symmetric matrix (N x N).
 *   - int N: Dimension of the matrices.
 *
 * OUTPUT:
 *   - None. The matrix W is updated in place with the normalized values.
 */
void get_norm(double **W, double **D, double **A, int N){
    int i,j;

    /*compute D' = D^(-1/2)
    ** D changes in place to D^(-1/2) */
    for (i=0;i<N;i++) {
        D[i][i] = 1/sqrt(D[i][i]);
    }

    /*Compute W=[D^(-1/2)]*A*[D^(-1/2)] :
    ** mark D' = D^(-1/2)
    since D' is diagonal we get: W[i][j] = A[i][j]*D'[i][i]*D'[j][j] */
    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            W[i][j]+=A[i][j]*D[i][i]*D[j][j];
        }
    }
}


/*
 * FUNCTION: matrix_multiply
 * -------------------------
 * Multiplies two matrices A and B and updates the resulting matrix C in-place.
 *
 * INPUTS:
 *   - double **C: Pointer to the result matrix (m x n)
 *   - double **A: Pointer to the first matrix (m x k).
 *   - double **B: Pointer to the second matrix (k x n).
 *   - int m: Number of rows in matrix A.
 *   - int k: Number of columns in matrix A and rows in matrix B.
 *   - int n: Number of columns in matrix B.
 */
void matrix_multiply(double **C,double **A, double **B, int m, int k, int n) {
    int i, j, p;
    
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            C[i][j] = 0.0;
            for (p = 0; p < k; p++) {
                C[i][j] += A[i][p] * B[p][j];
            }
        }
    }
}


/*
 * FUNCTION: matrix_transpose
 * --------------------------
 * Computes the transpose of a given matrix A and updates the resulting matrix T in-place.
 *
 * INPUTS:
 *   - double **T: Pointer to the result matrix (cols x rows)
 *   - double **A: Pointer to the original matrix (rows x cols).
 *   - int rows: Number of rows in matrix A.
 *   - int cols: Number of columns in matrix A.
 *
 */
void matrix_transpose(double **T,double **A, int rows, int cols) {
    int i, j;
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            T[j][i] = A[i][j];
        }
    }
}


/*
 * FUNCTION: Frobenius_norm
 * ------------------------
 * Computes the squared Frobenius norm of the difference between two matrices A and B.
 *
 * INPUTS:
 *   - double **A: Pointer to the first matrix (rows x cols).
 *   - double **B: Pointer to the second matrix (rows x cols).
 *   - int rows: Number of rows in the matrices.
 *   - int cols: Number of columns in the matrices.
 *
 * OUTPUT:
 *   - Returns the squared Frobenius norm of (A - B).
 */
double Frobenius_norm(double** A, double** B, int rows, int cols) {
    double sum = 0.0;
    int i, j;

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            double diff = A[i][j] - B[i][j];
            sum += diff * diff;
        }
    }
    return sum;
}


/*
 * FUNCTION: copy_matrix
 * ---------------------
 * Copies the contents of one matrix to another.
 *
 * INPUTS:
 *   - double **source: Pointer to the source matrix (n x k).
 *   - double **destination: Pointer to the destination matrix (n x k), modified in place.
 *   - int n: Number of rows in the matrices.
 *   - int k: Number of columns in the matrices.
 *
 * OUTPUT:
 *   - None. The destination matrix is updated with values from the source matrix.
 */
void copy_matrix(double** source, double** destination, int n, int k) {
    int i,j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            destination[i][j] = source[i][j];
        }
    }
}


/*
 * FUNCTION: make_step
 * -------------------
 * Performs one update step in the optimization of matrix H.
 * Updates new_H with the next iteration values based on the multiplicative update rule.
 *
 * INPUTS:
 *   - double **new_H: Pointer to the matrix where the updated values of H will be stored (n x k), modified in place.
 *   - double **W: Pointer to the input matrix W (n x n).
 *   - double **H_just_before_the_step: Pointer to the current estimate of H before the update (n x k).
 *   - int n: Number of rows in H and W.
 *   - int k: Number of columns in H.
 *
 * OUTPUT:
 *   - Error flag: 0 indicates success, 1 indicates failure
 */

int make_step(double** new_H ,double** W,double** H_just_before_the_step, int n, int k){
    int i,j;
    double num, denom;
    double** numerator, ** Htranspose, ** HtransposeH, ** denominator; 
    numerator=create_mat(n,k);
    if (numerator==NULL) {
        return 1;}
    Htranspose=create_mat(k,n);
    if (Htranspose==NULL) {
        free_2D_array(numerator,n);
        return 1;}
    HtransposeH=create_mat(k,k); 
    if (HtransposeH==NULL) {
        free_2D_array(numerator,n);
        free_2D_array(Htranspose,k);
        return 1;}
    denominator=create_mat(n,k);
    if (denominator==NULL) {
        free_2D_array(numerator,n);
        free_2D_array(Htranspose,k);
        free_2D_array(HtransposeH,k);
        return 1;}
    matrix_multiply(numerator,W, H_just_before_the_step, n, n, k);
    matrix_transpose(Htranspose,H_just_before_the_step, n, k);
    matrix_multiply(HtransposeH,Htranspose, H_just_before_the_step, k, n, k);
    matrix_multiply(denominator,H_just_before_the_step, HtransposeH, n, k, k );
    free_2D_array(Htranspose, k);
    free_2D_array(HtransposeH, k);
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            num = numerator[i][j];
            denom = denominator[i][j];
            if (denom == 0.0) {
                denom += 0.000001;
            }          
            new_H[i][j] = H_just_before_the_step[i][j] * (0.5 + 0.5 * (num / denom));}}
    free_2D_array(numerator, n);
    free_2D_array(denominator, n);   
    return 0;
}


/*
 * FUNCTION: get_symnmf
 * --------------------
 * Performs Symmetric Non-negative Matrix Factorization (SymNMF) on matrix W.
 * Updates new_H with the result of the SymNMF optimization.
 *
 * INPUTS:
 *   - double **new_H: Pointer to the matrix to store the optimized H (n x k), modified in place.
 *   - double **W: Pointer to the input symmetric matrix (n x n).
 *   - double **H: Pointer to the initial H matrix (n x k).
 *   - int n: Number of rows in H and W.
 *   - int k: Number of columns in H.
 * OUTPUT:
 *   - Error flag: 0 indicates success, 1 indicates failure
 */
int get_symnmf(double** new_H, double **W, double **H, int n, int k)
{
    int i,error_flag;
    double** H_just_before_the_step = create_mat(n,k);

    if (H_just_before_the_step==NULL) {
        return 1;
    }
    copy_matrix(H, H_just_before_the_step, n, k);
    
    for (i=0;i<MAX_ITER;i++)
    {
        error_flag=make_step(new_H, W, H_just_before_the_step, n, k);
            if (error_flag==1) {
                free_2D_array(H_just_before_the_step,n);
                return 1;
            }

        if (Frobenius_norm(new_H, H_just_before_the_step, n, k) < EPSILON)
        {
            free_2D_array(H_just_before_the_step, n);
            return 0;
        }
        copy_matrix(new_H, H_just_before_the_step, n, k);
    }
    free_2D_array(H_just_before_the_step, n);
    return 0;
}


