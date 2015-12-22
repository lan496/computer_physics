import numpy as np

def determinant(A):
    N = A.shape[0]
    sgn = 0

    for i in np.arange(N):
        max_elm_col = max([(abs(A[i][j]), j) for j in np.arange(i,N)])[1]
        for j in np.arange(N):
            A[j][i], A[j][max_elm_col] = A[j][max_elm_col], A[j][i]
        if max_elm_col != i:
            sgn += 1

        if A[i][i] == 0:
            return 0.0

        for j in np.arange(i+1,N):
            A[j] -= A[i] * (A[j][i] / A[i][i])

    det = (-1)**sgn * reduce(lambda x,y: x*y, [A[i][i] for i in np.arange(N)])
    return det

def main():
    A = np.array([[1,0,0],[0,1,2],[1,0,4]], dtype = float)

    print determinant(A)

main()


