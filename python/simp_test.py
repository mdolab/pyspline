from numpy import *
import sys
from matplotlib.pylab import plot,show
from pyspline import bintk,bvalu
x = linspace(0,10,12)
y = sin(0.5*x)

plot(x,y,'ko-')

k = 4
t = zeros(len(x)+k)

t[0:3] = 0
t[-1] = 10
t[-2] = 10
t[-3] = 10

#we have n-k values left

t_temp = linspace(0,10,10)
t[3:3+10] = t_temp

bcoef = bintk(x,y,t,4)

for i in xrange(12):
    print 'y,b:',y[i],bcoef[i]

plot(x,y,'ko-')

bcoef[-2] = -.5
bcoef[-1] = 0 


plot(x,bcoef,'ks')


#reinterpolate

x2 = linspace(0,10,50)
y2 = zeros(50)

for i in xrange(50):
    y2[i],inbv = bvalu(t,bcoef,k,0,x2[i])

plot(x2,y2,'-')

show()
