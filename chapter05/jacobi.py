import numpy as np

def max_nondiagonal(A):
    N = A.shape[0]
    tmp = [0.0, (0,0)]
    for i in np.arange(N):
        for j in np.arange(i+1,N):
            if abs(A[i][j]) > tmp[0]:
                tmp = [abs(A[i][j]), (i,j)]
    return tmp

def jacobi(A):
    N = A.shape[0]
    U = np.identity(N)

    while(True):
        non_dgn = max_nondiagonal(A)
        if non_dgn[0] < 1e-11:
            break
        p, q = non_dgn[1]
        
        z = np.sqrt(4.0*(A[p][q]**2) + (A[p][p] - A[q][q])**2)
        c2 = abs(A[p][p] - A[q][q])/z
        cs = np.sqrt((1.0 + c2)/2.0)
        sn = np.sqrt((1.0 - c2)/2.0)
        if -A[p][q]*(A[p][p] - A[q][q]) < 0:
            sn *= -1

        A_tmp = np.array(A)
        for i in np.arange(N):
            A_tmp[i][p] = cs*A[i][p] - sn*A[i][q]
            A_tmp[i][q] = sn*A[i][p] + cs*A[i][q]
        for j in np.arange(N):
            A_tmp[p][j] = cs*A[p][j] - sn*A[q][j]
            A_tmp[q][j] = sn*A[p][j] + cs*A[q][j]
        A_tmp[p][q] = sn*cs*(A[p][p] - A[q][q]) + c2*A[p][q]
        A_tmp[q][p] = A_tmp[p][q]
        A_tmp[p][p] = (cs**2)*A[p][p] + (sn**2)*A[q][q] - 2*cs*sn*A[p][q]
        A_tmp[q][q] = (sn**2)*A[p][p] + (cs**2)*A[q][q] + 2*cs*sn*A[p][q]

        U_tmp = np.array(U)
        for i in np.arange(N):
            U_tmp[i][p] = U[i][p] * cs - U[i][q] * sn
            U_tmp[i][q] = U[i][p] * sn + U[i][q] * cs
        
        A = np.array(A_tmp)
        U = np.array(U_tmp)
   
    eigenvalue = [(A[i][i], i) for i in np.arange(N)]
    eigenvalue.sort()
    pivot = [e[1] for e in eigenvalue]
    
    A_res = np.identity(N)
    U_res = np.identity(N)
    for i in np.arange(N):
        for j in np.arange(N):
            A_res[i][j] = A[pivot[i]][pivot[j]]
            U_res[i][j] = U[pivot[i]][pivot[j]]
    
    return A_res, U_res

def test():
    A = np.array([[5,4,1,1], [4,5,1,1], [1,1,4,2], [1,1,2,4]], dtype = float)
    AA, U = jacobi(A)
    print AA
    print U
    print A * U - AA * U
test()
