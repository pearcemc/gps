from pylab import *

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

no = Xz.shape[0]
ns = Xf.shape[1]

Kzz = evalK(Xz,Xz,K)
Kzf = evalK(Xz,Xf,K)
Kfz = Kzf.T
Kff = evalK(Xf,Xf,K)

Ef = evalM(Xf,M) + Kfz * Kzz.I * ( z - evalM(Xz,M) )
Ef
