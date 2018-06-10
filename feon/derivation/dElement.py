# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------
from __future__ import division
import numpy as np
from .lagrange import *
from .base import dElementBase
from ..tools import gl_quad1d
from .integration import tri_quad,Cubtri,Triex,Hammer,tetra_quad,Zctetra
import time

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

class Line(dElementBase):
    def __init__(self,ntype = "lagrange"):
        dElementBase.__init__(self,ntype)
        self._dim = 1

    def calc_ncoords(self):
        self._ncoords = np.array([0]+[np.linalg.norm(self.coords[0]-c) for c in self.coords[1:]])
        self._order = self.ncoords.shape[0]
        
    def calc_nbase(self,nodes):
        self._nodes = nodes
        assert self.nop >= 2,"At least two nodes needed"
        self._coords = np.array([nd.coord for nd in self.nodes])
        self.calc_ncoords()
        self._volume = self.ncoords[1]
        if self.ntype is "lagrange":
            for i in range(self.nop):
                self._n_base.append(lagrange1d(self.ncoords,i))
        else:
            raise AttributeError("Unknwn shape function type %r"(self.ntype,))

    def calc_ke(self,**kwargs):
        self._ke = gl_quad1d(self.func,self.nop,(0,self.volume))


    def func(self,x):
        shape = self.L["shape"]
        B = np.zeros((shape[0],shape[1]*self.nop))
        for i in range(self.nop):
            for val in self.L["mapping"].items():
                j = i*shape[1]+val[0][1]
                k = val[0][0]
                B[k,j] = self.n_base[i](x,val[1])
        return np.dot(np.dot(B.T,self.D),B)
    
    def set_D(self,D):
        self._D = D
        
    def set_L(self,**kwargs):
        self._L["shape"] = kwargs["shape"]
        self._L["mapping"] = kwargs["mapping"]

class Triangle(dElementBase):
    def __init__(self,ntype = "lagrange" ):
        dElementBase.__init__(self,ntype)
        self._ac = None


    @property
    def ac(self):
        return self._ac

        
    def calc_ncoords(self):
        order = self.order
        one = [((order -i)/order,1-(order-i)/order,0) for i in range(order)]
        two = [(0,(order -i)/order,1-(order-i)/order) for i in range(order)]
        three = [(i/order,0,1-(i/order)) for i in range(order)]
        self._ncoords = np.vstack((one,two,three))

            
    
    def calc_nbase(self,nodes,order = 1):
        self._order = order
        self._nodes = nodes
        if self.order == 1:
            assert self.nop == 3
        elif self.order == 2:
            assert self.nop == 6
        elif self.order == 3:
            assert self.nop == 10
        else:
            raise AttributeError("Only there order supported by far")
        
        self._coords = np.array([nd.coord for nd in self.nodes])
        self._volume = area_of_tri(self.coords[0],self.coords[1],self.coords[2])
        if self.ntype is "lagrange":
            for i in range(self.nop):
                self._n_base.append(lagrange2d_tri(self.coords[:3],self.volume,self.order,i))

    def func(self,x):
        shape = self.L["shape"]
        B = np.zeros((shape[0],shape[1]*self.nop))
        for i in range(self.nop):
            for val in self.L["mapping"].items():
                j = i*shape[1]+val[0][1]
                k = val[0][0]
                B[k,j] = self.n_base[i]((x[0],x[1]),val[1])
        return np.dot(np.dot(B.T,self.D),B)

    def set_D(self,D):
        self._D = D
        
    def set_L(self,**kwargs):
        self._L["shape"] = kwargs["shape"]
        self._L["mapping"] = kwargs["mapping"]

    def calc_ke(self,method = "Hammer",n = 5,**kwargs):
        meth = method
        self._ke = tri_quad(self.func,self.coords[:3],method=eval(meth)(n))

class Tetrahegon(dElementBase):
    def __init__(self,ntype = "lagrange" ):
        dElementBase.__init__(self,ntype)

            
    def calc_nbase(self,nodes,order = 1):
        self._order = order
        self._nodes = nodes
        if self.order == 1:
            assert self.nop == 4
        elif self.order == 2:
            assert self.nop == 10
        else:
            raise AttributeError("Only two order supported by far")
        
        self._coords = np.array([nd.coord for nd in self.nodes])
        self._volume = volume_of_tetra(self.coords[0],self.coords[1],self.coords[2],self.coords[3])
        if self.ntype is "lagrange":
            for i in range(self.nop):
                self._n_base.append(lagrange3d_tetra(self.coords[:4],self.volume,self.order,i))

    def func(self,x):
        shape = self.L["shape"]
        B = np.zeros((shape[0],shape[1]*self.nop))
        for i in range(self.nop):
            for val in self.L["mapping"].items():
                j = i*shape[1]+val[0][1]
                k = val[0][0]
                B[k,j] = self.n_base[i]((x[0],x[1],x[2]),val[1])
        return np.dot(np.dot(B.T,self.D),B)

    def set_D(self,D):
        self._D = D
        
    def set_L(self,**kwargs):
        self._L["shape"] = kwargs["shape"]
        self._L["mapping"] = kwargs["mapping"]

    def calc_ke(self,method = "Zctetra",**kwargs):
        meth = method
        self._ke = tetra_quad(self.func,self.coords[:4],method=eval(meth)())
        
if __name__ =="__main__":
    from feon.sa import *
    nds  = [Node([0,0,0]),Node([0,-1.,0.]),Node([1.,-1.,0]),Node([0,-1.,1.])]
    e = Tetrahegon()
    e.calc_nbase(nds)
    E = 210e6
    nu = 0.3
    a = E/(1+nu)/(1-2*nu)
    b = 1.-nu
    c = (1.-2*nu)/2.
    D = a*np.array([[b,nu,nu,0.,0.,0.],
                    [nu,b,nu,0.,0.,0.],
                    [nu,nu,b,0.,0.,0.],
                    [0.,0.,0.,c,0.,0.],
                    [0.,0.,0.,0.,c,0.],
                    [0.,0.,0.,0.,0.,c]])
    e.set_L(shape = (6,3),mapping = {(0,0):(1,0,0),
                                     (1,1):(0,1,0),
                                     (2,2):(0,0,1),
                                     (3,0):(0,1,0),
                                     (3,1):(1,0,0),
                                     (4,1):(0,0,1),
                                     (4,2):(0,1,0),
                                     (5,0):(0,0,1),
                                     (5,2):(1,0,0),
                                     })
    e.set_D(D=D)
    e.calc_ke()
    
    
####    pass
##    nds = (Node(0,0),Node(1,0),Node(1,2))
####    nds1 = (Node(0,0),Node(0.5,0),Node(0.5,0.25),Node(0.25,0),Node(0.5,0.125),Node(0.25,0.125))
########    p2 = [[0,0],[1,0],[1,1],[0,1]]
######    e = Line()
######    e2 = Rect()
##    e.calc_nbase(nds)
##    e2.calc_nbase(p2)
##    e.set_L(shape = (1,1),mapping = {(0,0):1})
##    e.set_D(D=1)
##    e.calc_ke()
    
##    e2.set_L(shape = (3,2),index = {(0,0):(1,0),(1,1):(0,1),(2,0):(0,1),(2,1):(1,0)})
##    e2.set_D(1)
##    e2.func((0.5,0.5))
    
##    E = 210e6
##    nu = 0.3
##    t = 0.025
##    A = t*E/(1-nu**2)
##
##    start = time.clock()
##    e = Triangle()
##    e.calc_nbase(nds,order=1)
##    e.set_L(shape = (3,2),mapping = {(0,0):(1,0),
##                                     (1,1):(0,1),
##                                     (2,0):(0,1),
##                                     (2,1):(1,0)})
##    e.set_D(D=A*np.array([[1,nu,0],
##                          [nu,1,0],
##                          [0,0,(1-nu)/2.]]))
##    e.calc_ke(method = "Triex")
##    end = time.clock()
##    print end-start
##    from feon.sa import *
##    
##    start = time.clock()
##    e2 = Tri2D11S(nds,E,nu,t)
##    e2.calc_Ke()
##    end = time.clock()
##    print end -start
    
