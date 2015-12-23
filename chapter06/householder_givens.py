import numpy as np

def strum(A):
    N = A.shape[0]
    a = [A[i][i] for i in np.arange(N)]
    b = [A[i][i+1] for i in np.arange(N-1)]
    
    def L(x):
        f = [1.0, a[0] - x]
        for i in np.arange(2,N+1):
            f.append((a[i-1] - x)*f[i-1] - f[i-2]*b[i-2]**2)
        cnt = 0
        for i in np.arange(N):
            if f[i] * f[i+1] < 0 or (f[i] == 0 and i > 0 and f[i-1] * f[i+1] < 0):
                cnt += 1
        return cnt

    nrm = sum([abs(e) for e in a]) + 2.0 * sum([abs(e) for e in b])
    
    eigenvalues = []
    for i in np.arange(N):
        ub = nrm
        lb = -nrm
        while(ub - lb > 1e-8):
            mid = (ub + lb) / 2.0
            if L(mid) <= i:
                lb = mid
            else:
                ub = mid
        eigenvalues.append(ub)

    eigen = []
    for l in eigenvalues:
        u = np.ones(N)
        u[1] = (l - a[0]) / b[0]
        for i in np.arange(2,N):
            u[i] = ((l - a[i-1])*u[i-1] - b[i-2]*u[i-2])/b[i-1]
        u /= np.linalg.norm(u)
        eigen.append([l, u])

    return eigen

def householder(A):
    N = A.shape[0]
    P = np.identity(N)
    
    for k in np.arange(N-2):
        s = A[k][k+1] + np.sign(A[k][k+1]) * np.sqrt(sum([A[k][m]**2 for m in np.arange(k+1,N)]))
        p = np.zeros(N)
        p[k+1] = 1.0
        for m in np.arange(k+2,N):
            p[m] = A[k][m]/s
        alpha = 2.0/sum([e**2 for e in p])

        Pk = np.identity(N) - alpha * np.array([[p[i] * p[j] for j in np.arange(N)] for i in np.arange(N)])
        
        Ak = np.copy(A)
        Ak[k][k] = A[k][k] * 1.0
        Ak[k][k+1] = -np.sign(A[k][k+1]) * np.sqrt(sum([A[k][m]**2 for m in np.arange(k+1,N)]))
        for m in np.arange(k+2,N):
            Ak[m][k] = Ak[k][m] = 0.0
        Ak[k+1][k] = Ak[k][k+1] * 1.0
        Ak[k+1:,k+1:] = np.array(Pk[k+1:,k+1:].dot(A[k+1:,k+1:].dot(Pk[k+1:,k+1:])))

        A = np.array(Ak)
        P = np.array(P.dot(Pk))

    return (A, P)

def householder_givens(A):
    A3, P = householder(A)
    eigen = strum(A3)
    for e in eigen:
        e[1] = np.dot(P, e[1])
    return eigen

def test():
    A = np.array([[5,4,1,1], [4,5,1,1], [1,1,4,2], [1,1,2,4]], dtype=float)
    eigen = householder_givens(A)
    for e in eigen:
        print e

test()
