from pylab import *
from copy import deepcopy

def meanf(Xi,c1,c2):
  if (-1)**Xi[0,1] <0:
    return c1
  return c2

def makeM(c1,c2):
  return lambda i: meanf(i,c1,c2)

def evalM(A,M):
  d = A.shape[0] 
  B = matrix(zeros((d,1)))
  for i in range(d):
    B[i] = M(A[i,:])
  return B 

def kernf(Xi,Xj,g,s1,s2,s3):
  if (-1)**(Xi[0,1] + Xj[0,1])>0:
    dl = 0
    if Xi[0,0]==Xj[0,0]:
      dl = s2
    return s3*exp(-g*(Xi[0,0] - Xj[0,0])**2) + s1 + dl 
  return s1

def makeK(g,s1,s2,s3):
  return lambda i,j: kernf(i,j,g,s1,s2,s3)

def evalK(A,B,K):
  d1 = A.shape[0]
  d2 = B.shape[0]
  C = matrix(zeros((d1,d2)))
  for i in range(d1):
    for j in range(d2):
      C[i,j] = K(A[i,:],B[j,:])
  return C

def plot_xy(E,color="blue"):
  B = deepcopy(E)
  if not B.shape[1] == 2:
    B.shape = (E.shape[0]/2,2)
  plot( array(B[:,0].T)[0],array(B[:,1].T)[0], color=color )

c1 = 5.
c2 = 5.
g = .05
s1 = 1. 
s2 = 1.
s3 = 5.

M = makeM(c1,c2)
K = makeK(g,s1,s2,s3)

z = matrix("1.;1.;1.01;1.01;1.05;1.05")
Xz = matrix("0.1,0;0.1,1.;0.9,0;0.9,1;1.3,0;1.3,1")

N = 20
T = linspace(0.10001,5,N)
Xf = matrix(zeros((2*N,2)))
for i in range(N):
  Xf[2*i,:]= matrix([T[i],0])
  Xf[2*i + 1,:]= matrix([T[i],1])


Kzz = evalK(Xz,Xz,K)
Kzf = evalK(Xz,Xf,K)
Kfz = Kzf.T
Kff = evalK(Xf,Xf,K)

Ef = evalM(Xf,M) + Kfz * Kzz.I * ( z - evalM(Xz,M) )
Ef

Vf = Kff - Kfz * Kzz.I * Kzf

qchisq = 2.448
Cn = 100
C = matrix(zeros((Cn,2)))
T = linspace(0,2*pi,Cn)
for i in range(Cn):
  C[i,0] = cos(T[i])
  C[i,1] = sin(T[i])

def do_curve(V,i):
  eva, eve = eig(Vf[2*i:2*i + 2,2*i:2*i + 2])
  return C * eve * diag(eva**.5) * qchisq

def make_stack(v,n):
  R = matrix(zeros((n,2)))
  for i in range(n):
    R[i,:] = v
  return R

plot_xy(Ef)
for i in range(Xf.shape[0]/2):
  L = do_curve(Vf,i) + make_stack(Ef[2*i:2*i + 2,:].T,Cn)
  plot_xy(L,color='red') 
