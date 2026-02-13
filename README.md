# inverse-compound
An algorithm to recover the original matrix from its kth-compound. We present a MATLAB implementation of our algorithm.
Suppose for an unknown real matrix $A$ with known dimension, we have access to its $k^{th}$ compound matrix $C_k(A)$ for a known $k$. We present an algorithm, which recovers $A$ from a given $C_k(A)$. 

We demonstrate two approaches for the recovery:

1. Eigen decomposition: Run the program **reconstruct_eigendecomp.m** to recover a square matrix which is diagonalizable. It is based on the compound operation on the similarity relation $A=U\Lambda U^{-1}$, where $\Lambda$ is a diagonal matrix with the eigenvalues on the diagonal and $U$ is an eigenvector matrix with columns as independent eigenvectors. Our initial version supports only real and simple eigenvalues, i.e., each eigenvalue has algebraic multiplicity one, and real eigenvectors. However, the latter version is planned to extend to the general condition.
   
2. SVD: Run the program **reconstruct_svd.m** to recover a real rectagular matrix. It operates on the SVD decomposition of the matrix $A=U\Sigma W^\top$, where $Sigma$ is a diagonal matrix with non-negative singular values on the diagonal. $U$ and $W$ are orthonormal and have left and right singular vectors in their columns, respectively.

# Set-up and Instructions
To test the program, run the above MATLAB program files respectively. Provide your Input matrix $A$ in the input section provided in both codes. Compare the recovered matrix $A_{rec}$ against the true matrix

