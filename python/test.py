from numpy import *
import sys
from matplotlib.pylab import plot,show
from pyspline import *

#Debugging script for testing derivative  constrained 2D b-spline surfaces
# N = 34
# x0 = 0.5*(1-cos(linspace(0,pi,N)))
# x0 = array([0.00000,0.00018,0.00255,0.00954,0.02090,0.03650,0.05640,0.08030,0.10801,0.13934,0.17395,0.21146,0.25149,0.29361,0.33736,0.38228,0.42820,0.47526,0.52324,0.57161,0.61980,0.66724,0.71333,0.75749,0.79915,0.83778,0.87287,0.90391,0.93072,0.95355,0.97251,0.98719,0.99668,1.00000])
# y0 = array([0.00000,0.00159,0.00748,0.01640,0.02600,0.03580,0.04560,0.05520,0.06430,0.07290,0.08070,0.08760,0.09340,0.09810,0.10133,0.10294,0.10249,0.10005,0.09610,0.09090,0.08490,0.07820,0.07100,0.06340,0.05570,0.04800,0.04030,0.03260,0.02480,0.01700,0.00982,0.00431,0.00103,0.00000])

# #Top surface of a naca0012 airfoil
# y0= 0.12/0.2*(0.2969*x0**0.5-0.1260*x0-0.3516*x0**2 + 0.2843*x0**3-0.1015*x0**4)


N = 161
x = linspace(0,1,N)
#x = array([0,0.0001,0.0325,0.05,0.1,0.2,0.35,0.5,0.65,0.85,0.925,1])
x = 0.5*(1-cos(linspace(0,pi,N)))

xu = x
xl = x


yu = 0.12/0.2*(0.2969*xu**0.5-0.1260*xu-0.3516*xu**2 + 0.2843*xu**3-0.1015*xu**4)
yl = -0.12/0.2*(0.2969*xl**0.5-0.1260*xl-0.3516*xl**2 + 0.2843*xl**3-0.1015*xl**4)

x = hstack([xu,xl[1:]])
y = hstack([yu,yl[1:]])
x = xu
y = yu

#create s
#N =2*N-1

s = zeros(N)

for i in xrange(N-1):
    s[i+1] = s[i] + sqrt((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)

s = s/s[-1] #Normalize s

iopt = -1
w = ones(N)
w[0] = 100
w[-1] = 100

k = 3
n0 = 15
smooth = 1e-5
nest = len(x)+k+1
print 'input:'
print 'smooth:',smooth

# do a uniform knot vector
t = zeros(nest)
t[0:k+1] = s[0]
t[n0-k-1:n0] = s[-1]
temp = 0.5*(1-cos(linspace(0,2*pi,n0-2*k)))
t[k:n0-k] = 0.5*(1-cos(linspace(0,pi,n0-2*k)))
#t[k:n0-k] = [0, 0.25,0.35,0.4,0.425,0.45,0.475, 0.5,0.525,0.55, 0.575,0.6,0.65,0.75,1]
#print 't',t 
#t[k:n0-k] = linspace(0,1,n0-2*k)

print 't_before:',t[0:n0]

n,tx,cx,fp,ier = curfit(iopt,s,x,w,s[0],s[-1],k,smooth,n0,t)
n,ty,cy,fp,ier = curfit(iopt,s,y,w,s[0],s[-1],k,smooth,n0,t)

cx[0] = 0
cy[0] = 0
cx[1] = 0
print 'n:',n
print 't:',t[0:n]
print 'cx:',cx[0:n-k-1]
print 'cy:',cy[0:n-k-1]

print 'fp:',fp
print 'ier:',ier


plot(x,y,'ko-')
plot(cx[0:n-k-1],cy[0:n-k-1],'s')


#Check interpolation
x_interp,ier = splev(tx[0:n],cx[0:n],k,s)
y_interp,ier = splev(ty[0:n],cy[0:n],k,s)


plot(x_interp,y_interp,'-')

show()
