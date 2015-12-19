import numpy as np

def gauss_seidel(A, b):
    N = b.size
    x_prev = np.array([b[i]/A[i][i]  for i in np.arange(N)])
    x_now = x_prev

    while(True):
        for i in np.arange(N):
            x_now[i] = b[i]
            for j in np.arange(i):
                x_now[i] -= A[i][j] * x_now[j]
            for j in np.arange(i+1,N):
                x_now[i] -= A[i][j] * x_prev[j]
            x_now[i] /= A[i][i]

        if np.linalg.norm(x_now - x_prev) < 1e-8:
            break

        x_prev = x_now

    return [e2 for e1 in x_prev for e2 in e1]

def main():
    A = np.array([[1,0,0],[0,1,2],[1,0,4]], dtype = float)
    b = np.array([[5], [12], [9]], dtype = float)
    
    print gauss_seidel(A, b)

main()
