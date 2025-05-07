import pandas as pd
import numpy as np
import sys
import math
from symnmf import SymnmfAlgorithm
from kmeans import KmeansAlgorithm
from sklearn.metrics import silhouette_score

def main():
    """
    Main function to perform clustering and evaluate results using SymNMF and KMeans.

    The function reads a dataset, extracts the number of clusters, performs clustering 
    using Symmetric Non-negative Matrix Factorization (SymNMF) and KMeans, and calculates 
    their silhouette scores.

    Arguments:
        None (Command-line arguments are used).

    Expected command-line arguments:
        - k: An integer representing the number of clusters (for SymNMF and KMeans).
        - data: A path to a CSV or TXT file containing the dataset.

    Returns:
        None
    """

    CASE_ERROR ="An Error Has Occurred\n"

    if len(sys.argv) != 3:
        print(CASE_ERROR)                 
        exit()

    data = sys.argv[2]
    if not (data.endswith((".txt", ".csv"))):
        print(CASE_ERROR)                 
        exit()
    points = pd.read_csv(data, header=None).values.tolist()
    N = len(points)             
    
    k = int(sys.argv[1])
    if (k<=1 or k>=N):
        print(CASE_ERROR)                 
        exit()

    symnmf_labels = SymnmfAlgorithm(points, k)
    if symnmf_labels==None:
        print(CASE_ERROR)  
        exit()

    symnmf_score = silhouette_score(points, symnmf_labels)
    Kmeans_labels = KmeansAlgorithm(points, k)
    Kmeans_score = silhouette_score(points, Kmeans_labels)
    print("nmf:","{:.4f}".format(symnmf_score))
    print("kmeans:","{:.4f}".format(Kmeans_score))



if __name__ == "__main__":
    main()
  