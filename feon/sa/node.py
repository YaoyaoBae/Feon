# -*- coding: utf-8 -*-
# ------------------------------------
#  Author: YAOYAO PEI
#  E-mail: yaoyao.bae@foxmail.com
#  License: Hubei University of Technology License
# -------------------------------------
from ..base import NodeBase
import numpy as np

#~ *********************************************************************
#~ ****************************  CLASS NODE ****************************
#~ *********************************************************************
class Node(NodeBase):

    def __init__(self,*coord):
        NodeBase.__init__(self,*coord)
        self._dof = len(self.nAk)
        self._disp = dict.fromkeys(self.nAk,0.)
        self._force = dict.fromkeys(self.nBk,0.)


    @property
    def force(self):
        return self._force
        
    @property
    def disp(self):
        return self._disp
        
    def init_keys(self):
        """
        K*D = F
        K*A = B
        """
        if self.dim == 2:
            self.set_nAk(("Ux","Uy","Phz"))
            self.set_nBk(("Fx","Fy","Mz"))
        elif self.dim == 3:
            self.set_nAk(("Ux","Uy","Uz","Phx","Phy","Phz"))
            self.set_nBk(("Fx","Fy","Fz","Mx","My","Mz"))

    def init_unknowns(self,*unknowns):
        for key in unknowns:
            if key in self.nAk:
                self._disp[key] = None
            else:
                raise AttributeError("Unknow disp name(%r)"%(unknowns,))
    
    def set_force(self,**forces):
        for key in forces.keys():
            if key in self.nBk:
                self._force[key] += forces[key]
            else:
                raise AttributeError("Unknow focre name(%r)"%(forces,))
    
    def clear_force(self):
        for key in self.nBk:
            self._force[key] = 0.

    def get_force(self):
        return self._force

    def set_disp(self,**disp):
        for key in disp.keys():
            if key in self.nAk:
                self._disp[key] = disp[key]
            else:
                raise AttributeError("Unknow disp name(%r)"%(disp,))
            
    def clear_disp(self):
        for key in self.nAk:
            self._disp[key] = 0.
        
    def get_disp(self):
        return self._disp
    
        
if __name__ == "__main__":
    n = Node(1,2)
