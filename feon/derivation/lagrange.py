# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------
from __future__ import division
import numpy as np
from functools import reduce
from collections import Counter
import mpmath as mp

def area_of_tri(p1,p2,p3):
    x,y = p1[0],p1[1]
    xj,yj = p2[0],p2[1]
    xm,ym = p3[0],p3[1]
    return 0.5*((xj*ym-yj*xm) +(yj-ym)*x +(xm-xj)*y)

def volume_of_tetra(p1,p2,p3,p4):
    V = np.ones((4,4))
    V[0,1:],V[1,1:],V[2,1:],V[3,1:]  = p1,p2,p3,p4
    return abs(np.linalg.det(V)/6.)

def differentiation1d(func):
    mp.mp.dps = 15;mp.mp.pretty = True
    def _derivate(x,deriv = 0):
        ret = mp.diff(func,x,deriv)
        return float(ret)
    return _derivate

def differentiation2d(func):
    mp.mp.dps = 15;mp.mp.pretty = True
    def _derivate(x,deriv = (0,0)):
        ret = mp.diff(func,(x[0],x[1]),deriv)
        return float(ret)
    return _derivate

def differentiation3d(func):
    mp.mp.dps = 15;mp.mp.pretty = True
    def _derivate(x,deriv = (0,0,0)):
        ret = mp.diff(func,(x[0],x[1],x[2]),deriv)
        return float(ret)
    return _derivate

def lagrange1d(coords,n):
    pv = coords[n]
    p = coords
    nx = p.shape[0]
    @differentiation1d
    def shape_func(x):
        l = ((x-p[i])/(pv-p[i]) for i in range(nx) if p[i] != pv)
        return reduce(lambda a,b:a*b,l)
    
    return shape_func

def lagrange2d_rect(coords,n):
    pi,pj = np.unique(coords[:,0]),np.unique(coords[:,1])
    nx,ny = pi.shape[0],pj.shape[0]
    pv = coords[n]
    
    @differentiation2d
    def shape_func(x,y):
        lx = ((x-pi[i])/(pv[0]-pi[i]) for i in range(nx) if pi[i] != pv[0])
        ly = ((y-pj[i])/(pv[1]-pj[i]) for i in range(nx) if pj[i] != pv[1])
        l = map(lambda a,b:a*b,zip(lx,ly))
        return reduce(lambda a,b:a*b,l)
    return shape_func

def lagrange2d_tri(coords,volume,order,n):
    mc = coords
    m = order
    v = volume
    n = n
    
    @differentiation2d
    def shape_func(x,y):
        L1,L2,L3 = area_of_tri(np.array([x,y]),mc[1],mc[2])/v,\
                   area_of_tri(np.array([x,y]),mc[2],mc[0])/v,\
                   area_of_tri(np.array([x,y]),mc[0],mc[1])/v  
        if m == 1:
            L = [L1,L2,L3]
            return L[n]
        if m == 2:
            L = [2*(L1-1)*L1,2*(L2-1)*L2,2*(L3-1)*L3,
                 4*L1*L2,4*L2*L3,4*L3*L1]
            return L[n]
        if m == 3:
            L = [1/2.*(3*L1-1)*(3*L1-2)*L1,1/2.*(3*L2-1)*(3*L2-2)*L2,1/2.*(3*L3-1)*(3*L3-2)*L3,
                 9/2.*L1*l2*(3*L1-1),9/2.*L1*l2*(3*L2-1),
                 9/2.*L2*l3*(3*L2-1),9/2.*L2*l3*(3*L3-1),
                 9/2.*L1*l3*(3*L3-1),9/2.*L1*l3*(3*L1-1),
                 27*L1*L2*L3]
            return L[n]
    return shape_func


def lagrange3d_tetra(coords,volume,order,n):
    mc = coords
    m = order
    v = volume
    n = n
    
    @differentiation3d
    def shape_func(x,y,z):
        L1,L2,L3,L4= volume_of_tetra(np.array([x,y,z]),mc[1],mc[2],mc[3])/v,\
                     volume_of_tetra(mc[0],np.array([x,y,z]),mc[2],mc[3])/v,\
                     volume_of_tetra(mc[0],mc[1],np.array([x,y,z]),mc[3])/v,\
                     volume_of_tetra(mc[0],mc[1],mc[2],np.array([x,y,z]))/v
        
        if m == 1:
            L = [L1,L2,L3,L4]
            return L[n]
        if m == 2:
            L = [2*(L1-1)*L1,2*(L2-1)*L2,2*(L3-1)*L3,2*(L4-1)*L4,
                 4*L1*L2,4*L1*L3,4*L1*L4,4*L2*L3,4*L3*L4,4*L2*L4]
            return L[n]
    return shape_func
    

if __name__ =="__main__":
    mc = np.array([[0,0,0],[0,-1.,0.],[1.,-1.,0],[0,-1.,1.]])
    def f(x,y,z):
        return volume_of_tetra(np.array([x,y,z]),mc[1],mc[2],mc[3])
    def fx(x,y,z):
        return x*y+z
    a = mp.diff(fx,(0.1,1,0.5),(0,1,0))

