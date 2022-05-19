#!/usr/bin/python

import numpy

xvals = (-1,1,1)
sigmavals = (0.1,0.1,0.1)
npar = 2

F = numpy.zeros([npar,npar])

x = zip(xvals,sigmavals)
print(*x)

for x,sigma in zip(xvals,sigmavals):
    for i in range(npar):
        if i==0:
            dfdpi = x
        else:
            dfdpi = 1
        for j in range(npar):
            if j==0:
                dfdpj = x
            else:
                dfdpj = 1
                
            F[i,j] += sigma**-2*dfdpi*dfdpj
    
print(numpy.mat(F).I) # invert the matrix
