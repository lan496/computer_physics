import numpy as np

def gaussian_elimation(A, b):
    N = b.size
    pivot = np.arange(N)
    B = np.hstack((A, b))

    for i in np.arange(N):
        max_elm_col = max([(abs(B[i][j]), j) for j in np.arange(i,N)])[1]
        for j in np.arange(N):
            B[j][i], B[j][max_elm_col] = B[j][max_elm_col], B[j][i]
        pivot[i], pivot[max_elm_col] = pivot[max_elm_col], pivot[i]

        for j in np.arange(i+1,N):
            B[j] -= B[i] * (B[j][i] / B[i][i])

    x = np.zeros(N)
    for i in np.arange(N-1,-1,-1):
        x[i] = B[i][N]
        for j in np.arange(i+1,N):
            x[i] -= B[i][j] * x[j]
        x[i] /= B[i][i]

    return x[pivot]

def main():
    A = np.array([[1,0,0],[0,1,2],[1,0,4]], dtype = float)
    b = np.array([[5], [12], [9]], dtype = float)

    print gaussian_elimation(A, b)

main()

