import pandas as pd
import numpy as np
import sys
import symnmfmodule as sm
np.random.seed(1234)



def create_initial_H(matrix, rows, cols):
    """
    Creates an initial H matrix for NMF.

    @param matrix Input matrix (used to compute the mean value)
    @param rows Number of rows in the resulting H matrix
    @param cols Number of columns in the resulting H matrix
    @return A 2D list representing the initial H matrix
    """
    matrix = np.asarray(matrix)
    mean_val = np.mean(matrix)
    limit = 2 * (mean_val / cols) ** 0.5
    H_matrix = np.random.uniform(0.0, limit, size=(rows, cols))
    return H_matrix.tolist()


def print_matrix(mat, row_count, col_count):
    """
    Prints a matrix with formatted values.

    Parameters:
        mat (list of list of float): The matrix to print.
        row_count (int): The number of rows in the matrix.
        col_count (int): The number of columns in the matrix.

    Returns:
        None
    """
    for row_idx in range(row_count):
        line = ",".join("{:.4f}".format(mat[row_idx][col_idx]) for col_idx in range(col_count))
        print(line)
        

def SymnmfAlgorithm(points, k):
    """
    Performs Symmetric Non-negative Matrix Factorization (SymNMF) and assigns cluster labels.

    Parameters:
        points (list of list of floats): The input data points (matrix).
        k (int): The number of clusters or components.

    Returns:
        list of int: A list of cluster labels (0 to k-1), where each label corresponds to the assigned cluster, OR None in case of failure

    """

    N = len(points)
    W = sm.norm(points)
    if W==None:
        return None        
    H0 = create_initial_H(W, N, k)
    H = sm.symnmf(W, H0, k)
    if H==None:
        return None
    labels = np.argmax(np.array(H), axis=1).tolist()
    return labels


def check_input(argv):
    """
    Validates command-line arguments and loads the dataset.

    The function checks that input format is correct, loads the dataset, and extracts the required parameters.

    Arguments:
        argv: List of command-line arguments.

    Returns:
        - flag_error (int): 1 if input is invalid, 0 otherwise.
        - k (int or None): Number of clusters.
        - goal (str or None): One of ['symnmf', 'sym', 'ddg', 'norm'].
        - points (list of list of float or None): Loaded dataset."""

    
    if len(sys.argv) != 4:                 
        return 1, None, None, None

    data = sys.argv[3]
    if not (data.endswith((".txt", ".csv"))):
        return 1, None, None, None

    goal = sys.argv[2]
    if goal not in ["symnmf","sym","ddg","norm"]:
        return 1, None, None, None
    
    data = sys.argv[3]
    points = pd.read_csv(data, header=None).values.tolist()
    N = len(points)
    
    k = int(sys.argv[1])
    if (goal == "symnmf") and (k<=1 or k>=N):
        return 1, None, None, None

    return 0, k, goal, points


def main():
    
    """
    Main function to validate input and run the selected matrix computation.

    Arguments:
        None (command-line arguments are used).

    Expected command-line arguments:
        - k: Integer, number of clusters (used only for 'symnmf').
        - goal: One of ['symnmf', 'sym', 'ddg', 'norm'].
        - data: Path to a .txt or .csv file containing the dataset.

    Returns:
        None
    """
    CASE_ERROR ="An Error Has Occurred\n"
    flag_error, k, goal, points = check_input(sys.argv)
    if flag_error == 1:
        print(CASE_ERROR)
        exit()
    
    N = len(points)
    # PROGRAM
    if (goal == "sym"):
        A = sm.sym(points)
        if A==None:
            print(CASE_ERROR)
            exit()
        print_matrix(A, N, N)
    if (goal == "ddg"):
        D = sm.ddg(points)
        if D==None:
            print(CASE_ERROR)
            exit()
        print_matrix(D, N, N)
    if (goal in ["norm", "symnmf"]):
        W = sm.norm(points)
        if W==None:
            print(CASE_ERROR)
            exit()
    if (goal == "symnmf"):
        H0 = create_initial_H(W, N, k)
        H = sm.symnmf(W, H0, k)
        if H==None:
            print(CASE_ERROR)
            exit()
        print_matrix(H, N, k)
    if (goal == "norm"):
        print_matrix(W, N, N)


if __name__ == "__main__":
    main()
  