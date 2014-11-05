from pylab import *
from copy import deepcopy

#def meanf(t,x,r,cx,cy,x0,y0):
#  return x*(cx - sqrt(cx**2 + x0**2 -2*x0*cx -2*r*t)) + (1-x)*(cy - sqrt(cy**2 + y0**2 -2*y0*cy -2*r*t))

meanf = lambda t,x,x0,y0,cx,cy,r: x*((x0 - cx)*exp(-r*t) + cx) + (1-x)*((y0 - cy)*exp(-r*t) + cy)

def makeM(x0,y0,cx,cy,r):
  return lambda t,x: meanf(t,x,x0,y0,cx,cy,r)

def evalM(A,M):
  d = A.shape[0] 
  B = matrix(zeros((d,1)))
  for i in range(d):
    B[i] = M(A[i,0],A[i,1])
  return B 

def kernf(ti,xi,tj,xj,g,s1,s2,s3):
  if (-1)**(xi + xj)>0:
    dl = 0
    if ti == tj:
      dl = s2
    return s3*exp(-g*(ti - tj)**2) + dl + s1
  return 0

def makeK(g,s1,s2,s3):
  return lambda ti,xi,tj,xj: kernf(ti,xi,tj,xj,g,s1,s2,s3)

def evalK(A,B,K):
  d1 = A.shape[0]
  d2 = B.shape[0]
  C = matrix(zeros((d1,d2)))
  for i in range(d1):
    for j in range(d2):
      C[i,j] = K(A[i,0],A[i,1],B[j,0],B[j,1])
  return C

def plot_xy(E,color="blue",scat=False):
  B = deepcopy(E)
  if not B.shape[1] == 2:
    B.shape = (E.shape[0]/2,2)
  if scat:
    scatter( array(B[:,1].T)[0],array(B[:,0].T)[0], color=color )
  else:
    plot( array(B[:,1].T)[0],array(B[:,0].T)[0], color=color )

def do_curve(V,i):
  eva, eve = eig(Vf[2*i:2*i + 2,2*i:2*i + 2])
  return C * eve * diag(eva**.5) * qchisq

def make_stack(v,n):
  R = matrix(zeros((n,2)))
  for i in range(n):
    R[i,:] = v
  return R

cx = 5.
cy = 4.
x0 = 1.
y0 = 2.
r = 0.4

g = .5
s1 = 0.05 
s2 = 0.01
s3 = 0.001

M = makeM(x0,y0,cx,cy,r)
K = makeK(g,s1,s2,s3)

#z = matrix([y0,x0,y0+.05,x0+.005,y0+.1,x0-.1,y0+.05,x0+.1,y0+.13,x0+.15]).T
#N = 5
#T = linspace(0.0,1.,N)
#Xz = matrix(zeros((2*N,2)))
#for i in range(N):
#  Xz[2*i,:]= matrix([T[i],0])
#  Xz[2*i + 1,:]= matrix([T[i],1])

z = matrix([[y0,x0]]).T
Xz = matrix("0,0;0,1")

N = 200
T = linspace(0.10001,10,N)
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
Tc = linspace(0,2*pi,Cn)
for i in range(Cn):
  C[i,0] = cos(Tc[i])
  C[i,1] = sin(Tc[i])

plot_xy(Ef)
for i in range(Xf.shape[0]/2):
  L = do_curve(Vf,i) + make_stack(Ef[2*i:2*i + 2,:].T,Cn)
  plot_xy(L,color='red') 


q = matrix(randn(Xf.shape[0],1))
eva, eve = eig(Vf)
Vfh = real(eve) * diag(real(eva)**.5)
Yf = Vfh * q + Ef

plot_xy(z,color="purple",scat=True)
plot_xy(Yf,color="green")
show()


