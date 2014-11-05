from pylab import *

def distance(v1,v2):
  return norm(v1 - v2)

def vhat(i,j,t,F,T,big=999):
  F = matrix(F)
  c = 0
  m = 999
  ij = matrix([i,j])
  nT = len(T)-1
  while T[c] < t:
    d = distance(F[c,:],ij)
    if d < m:
      m = d
    c += 1
    if c>nT:
      break
  return m
  
def image_meanf(i,j,t,h,l,F,T):
  return (1 + exp(h*vhat(i,j,t,F,T) + l))**(-1)

def makeL(h,l,F,T):
  return lambda i,j,t: image_meanf(i,j,t,h,l,F,T) 

h = 3.
l = 0.2
F = deepcopy(Yf)
F.shape = (Yf.shape[0]/2,2)
T = T

L = makeL(h,l,F,T)

xlb = min(F[:,0])[0,0] -3
xub = max(F[:,0])[0,0] +3
ylb = min(F[:,1])[0,0] -3
yub = max(F[:,1])[0,0] +3
Mp = 50
Ms = linspace(xlb,xub,Mp)
Np = 50
Ns = linspace(xlb,xub,Np)

#t = max(T)
t = 30.
I = matrix(zeros((Mp,Np)))
for i in range(Mp):
  for j in range(Np):
    I[i,j] = L(Ms[i],Ns[j],t)

fig, ax = plt.subplots(figsize=(6,6))
ax.imshow(I, interpolation='none', extent=[xlb,xub,yub,ylb], origin='lower')
