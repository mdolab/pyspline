from numpy import *
import sys
from matplotlib.pylab import plot,show
from pyspline import *
N = 15
x = linspace(0,1,N)
x = array([0,0.0001,0.0325,0.05,0.1,0.2,0.35,0.5,0.65,0.85,0.925,1])
x = 0.5*(1-cos(linspace(0,pi,N)))
xu = x[::-1]
xl = x



yu = 0.12/0.2*(0.2969*xu**0.5-0.1260*xu-0.3516*xu**2 + 0.2843*xu**3-0.1015*xu**4)
yl = -0.12/0.2*(0.2969*xl**0.5-0.1260*xl-0.3516*xl**2 + 0.2843*xl**3-0.1015*xl**4)

x = hstack([xu,xl[1:]])
y = hstack([yu,yl[1:]])

#create s 
s = zeros(2*N-1)

for i in xrange(2*N-1-1):
    s[i+1] = s[i] + sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)

s = s/s[-1] #Normalize s

# for i in xrange(2*N-1):
#     print 's,x,y:',s[i],x[i],y[i]

#plot(x,y,'k-')

k = 4
t = zeros(len(x)+k)

t[0:3] = 0
t[-1] = 1.0
t[-2] = 1.0
t[-3] = 1.0
t[-4] = 1.0
#we have n-k values left

t[4:4+(2*N-1)-4] = s[2:-2]

bcoef_x = bintk(s,x,t,4)
bcoef_y = bintk(s,y,t,4)
plot(bcoef_x,bcoef_y,'rs')


#Real Line
x2 = linspace(0,1,500)
x2 = 0.5*(1-cos(linspace(0,pi,500)))
y2 = 0.12/0.2*(0.2969*x2**0.5-0.1260*x2-0.3516*x2**2 + 0.2843*x2**3-0.1015*x2**4)
plot(x2,y2,'-')


#t,bcoef_x,n,k = bint4(s,x,1,2,0,0,1)
#t,bcoef_y,n,k = bint4(s,y,2,2,0,0,1)
#plot(bcoef_x,bcoef_y,'kd')
# print 't4:',t
# print 'bcoefx:',bcoef_x
# print 'bceofy:',bcoef_y
# print 'n:',n
# print 'k:',k


#Re interpolate
s = linspace(0,1,500)
x3 = zeros(500)
y3 = zeros(500)
for i  in xrange(len(x2)):
    x3[i],inbv = bvalu(t,bcoef_x,4,0,s[i])
    y3[i],inbv = bvalu(t,bcoef_y,4,0,s[i])

plot(x3,y3)

show()
