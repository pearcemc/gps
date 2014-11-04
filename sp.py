from pylab import *

def meanf(i,j,r,s,g,l):
  return 1./(g*sqrt((i-r)**2 + (s-j)**2) + l)

#def makeM(r,s,g,l):
#  return lambda i,j: meanf(i,j,r,s,g,l)

def logif(i,j,r,s,g):
  return (1. + exp(-g*sqrt((i-r)**2 + (s-j)**2)))**(-1)

def makeM(r,s,g):
  return lambda i,j: logif(i,j,r,s,g)

def kernf(i,j,r,s,s1,s2,b):
  ep = 0
  if i==r and s==j:
    ep = s1
  return s2*exp(sqrt((i-r)**2 + (s-j)**2)/b) + ep

def makeK(s1,s2,b):
  return lambda i,j,r,s: kernf(i,j,r,s,s1,s2,b)

M = makeM(2.,2.,-1.)
K = makeK(0.1,1.,2.)

R = 50
S = linspace(0,50,R)

EY = matrix(zeros((R,R)))
for i in range(R):
  for j in range(R):
    EY[i,j] = M(S[i],S[j])

imshow(EY)
