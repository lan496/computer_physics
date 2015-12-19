import numpy as np

def inverse_matrix(A):
    N = A.ndim + 1
    B = np.hstack((A, np.identity(N)))

    for i in np.arange(N):
        for j in np.arange(i+1,N):
            B[j] -= B[i] * (B[j][i] / B[i][i])

    for i in np.arange(N-1,-1,-1):
        B[i] /= B[i][i]
        for j in np.arange(i):
            B[j] -= B[i] * B[j][i]

    print B
    inv = np.hsplit(B, 2)[1]
    return inv

def main():
    A = np.array([[1,0,0],[0,1,2],[1,0,4]], dtype = float)
    b = np.array([[5], [12], [9]], dtype = float)

    print inverse_matrix(A)

main()


