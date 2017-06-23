# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------
from __future__ import division
import numpy as np

def pair_wise(L,end = False):
    pool = tuple(L)
    n = len(pool)
    assert n >=2,"Length of iterable must greater than 3"           
    if end is True:
        indices = list(range(n))
        for i in indices:
            if i < n -1:
                yield (pool[i],pool[i+1])
            else:
                yield (pool[i],pool[0])
    if end is False:
        indices = list(range(n-1))
        for i in indices:
            yield (pool[i],pool[i+1])
            
def gl_quad1d(fun,n,x_lim = None,args =()):
    if x_lim is None:
        a,b = -1,1
    else:
        a,b = x_lim[0],x_lim[1]
    
    if not callable(fun):
        return (b-a)*fun
    
    else:
        loc,w = np.polynomial.legendre.leggauss(n)
        s = (1/2.*(b-a)*fun((b-a)*v/2.+(a+b)/2.,*args)*w[i]
             for i,v in enumerate(loc))
        return sum(s)

            
def gl_quad2d(fun,n,x_lim = None,y_lim = None,args=()):
    
    if x_lim is None:
        a,b= -1,1
    else:
        a,b= x_lim[0],x_lim[1]
    if y_lim is None:
        c,d = -1,1
    else:
        c,d  = y_lim[0],y_lim[1]

    if not callable(fun):
        return (b-a)*(d-c)*fun
    else:
        loc,w = np.polynomial.legendre.leggauss(n)
        s = (1/4.*(b-a)*(d-c)*fun(((b-a)*v1/2.+(a+b)/2.,
                               (d-c)*v2/2.+(c+d)/2.),*args)*w[i]*w[j]
             for i,v1 in enumerate(loc)
             for j,v2 in enumerate(loc))
        return sum(s)



def gl_quad3d(fun,n,x_lim = None,y_lim = None,z_lim = None,args=()):
    if x_lim is None:
        a,b = -1, 1
    else:
        a,b= x_lim[0],x_lim[1]
    
    if y_lim is None:
        c ,d = -1,1
    else:
        c ,d = y_lim[0],y_lim[1]
    if z_lim is None:
        e,f= -1,1
    else:
        e ,f = z_lim[0],z_lim[1]

    if not callable(fun):
        return (b-a)*(d-c)*(f-e)*fun
    else:
        loc,w = np.polynomial.legendre.leggauss(n)
        s = (1/8.*(b-a)*(d-c)*(f-e)*fun(((b-a)*v1/2.+(a+b)/2.,
                                     (d-c)*v2/2.+(c+d)/2.,
                                     (f-e)*v3/2.+(e+f)/2.),*args)*w[i]*w[j]*w[k]
             for i,v1 in enumerate(loc)
             for j,v2 in enumerate(loc)
             for k,v3 in enumerate(loc))
        return sum(s)

def fun2(x,a,b):
    return a*x[0]*x[1]*np.e**(b*x[2])
if __name__ == "__main__":
    pass


    

    
##    
