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
        eigen.append((l, u))

    return eigen

def test():
    A = np.array([[2,-1,0,0], [-1,2,-1,0], [0,-1,2,-1], [0,0,-1,2]])
    print strum(A)

test()
