1# -*- coding: utf-8 -*-
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
        self._flowrate = dict.fromkeys(self.nBk,0.)#volumetric flow rate
        self._head = dict.fromkeys(self.nAk,None)#potential fluid head


    @property
    def head(self):
        return self._head
        
    @property
    def flowrate(self):
        return self._flowrate
        
    def init_keys(self):
        self.set_nAk(("H"))
        self.set_nBk(("FR"))
            
    
    def set_flowrate(self,**flowrate):
        self._flowrate["FR"] = val
    
    def clear_flowrate(self):
        for key in self.nBk:
            self._flowrate[key] = 0.

    def get_flowrate(self):
        return self._flowrate

    def set_head(self,val):
        self._head["H"] = val
            
    def clear_head(self):
        for key in self.nAk:
            self._head[key] = 0.
        
    def get_head(self):
        return self._head
    
        
if __name__ == "__main__":
    n = Node(1,2)
