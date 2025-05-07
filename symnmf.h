# ifndef SYMNMF_H_
# define SYMNMF_H_

/* Frees memory allocated for a 2D array*/
void free_2D_array(double **array, int rows);

/* Computes a symmetric matrix A based on pairwise distances between points */
void get_sym(double **points, double **A, int N, int d);

/* Computes the diagonal degree matrix D from a symmetric matrix A */
void get_diag(double **D, double **A, int N);

/* Computes the normalized matrix W using the symmetric matrix A and diagonal matrix D */
void get_norm(double **W, double **D, double **A, int N);

/* Performs Symmetric Non-negative Matrix Factorization (SymNMF) on matrix W */
int get_symnmf(double** newH, double **W, double **H, int n, int k);

/* Allocates and initializes a 2D matrix of doubles with all elements set to zero */
double **create_mat(int num_of_rows, int num_of_cols);

# endif
