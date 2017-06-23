# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

import numpy as np

class dElementBase(object):
    def __init__(self,ntype = "lagrange"):
        self._ntype = ntype
        self._nodes = []
        self._coords = None
        self._ncoords = None
        self._dim = 1
        self._volume = 0.
        self._n_base = []
        self._L = dict.fromkeys(("shape","index"),None)
        self._D = None
        self._ke = None
        self._order = 1
        
    def __repr__(self):
        return "%s Element(order=%s)"%(self.__class__.__name__,self.order)


    @property
    def ntype(self):
        return self._ntype
    
    @property
    def nodes(self):
        return self._nodes

    @property
    def order(self):
        return self._order
    
    @property
    def ke(self):
        return self._ke

    @property
    def D(self):
        return self._D

    @property
    def L(self):
        return self._L
    
    @property
    def n_base(self):
        return self._n_base
    
    @property
    def coords(self):
        return self._coords

    @property
    def ncoords(self):
        return self._ncoords

    @property
    def dim(self):
        return self._dim

    @property
    def nop(self):
        return len(self.nodes)

    @property
    def volume(self):
        return self._volume

        
    
if __name__ =="__main__":
    e = dElementBase("natrural")
