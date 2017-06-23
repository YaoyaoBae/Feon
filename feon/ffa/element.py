# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------

from __future__ import division
import numpy as np
from ..base import ElementBase
#~ *********************************************************************
#~ ************************* CLASS STRUCTELEMENT ***********************
#~ *********************************************************************


class Element(ElementBase):

    def __init__(self,nodes):
        ElementBase.__init__(self,nodes)
        self.init_nodes(nodes)
        self._velocity = dict.fromkeys(self.eIk,0.)
        
    def init_nodes(self,nodes):
        v = np.array(nodes[0].coord)-np.array(nodes[1].coord)
        le = np.linalg.norm(v)
        self._volume = le
        

    def init_keys(self):
        if self.dim == 2:
            self.set_eIk(("Vx","Vy"))
        if self.dim == 3:
            self.set_eIk(("Vx","Vy","Vz"))
    
    @property
    def T(self):
        return self._T

    
    @property
    def ke(self):
        return self._ke

    @property
    def velocity(self):
        return self._velocity

    
    def calc_ke(self):
        pass

    def calc_T(self):
        pass
        
    def calc_Ke(self):
        self.calc_T()
        self.calc_ke()
        self._Ke = np.dot(np.dot(self.T.T,self.ke),self.T)
        
        
    def evaluate(self):
        pass
        

    def distribute_velocity(self):
        n = len(self.eIk)
        for i,val in enumerate(self.eIk):
            self._velocity[val] += self._undealed_velocity[i::n]


class E1D(Element):

    def __init__(self,nodes,Kxx,A):
        Element.__init__(self,nodes)
        self.Kxx = Kxx
        self.A = A
        
        
    def init_keys(self):
        self.set_eIk(["Vx"])

    
    def calc_T(self):
        self._T = np.array([[1,0],
                            [0,1]])

    def calc_Ke(self):
        self.calc_T()
        self._Ke = _calc_ke_for_1d_e(self.Kxx,self.A,self.volume)


    def evaluate(self):
        h = np.array([[nd.head[key] for nd in self.nodes for key in nd.nAk[:self.ndof]]])
        a = np.array([[-1./self.volume,1./self.volume]])
        self._undealed_velocity = -self.Kxx*np.dot(a,h.T)
        self.distribute_velocity()
        
            
def _calc_ke_for_1d_e(Kxx,A,L):
    return np.array([[Kxx*A/L,-Kxx*A/L],
                     [-Kxx*A/L,Kxx*A/L]])

if __name__ == "__main__":
    pass


    
